"""Shared line-detection and empirical catalog-alignment statistics.

The circular-shift test evaluates whether a catalog's frequency pattern aligns
with detected spectral features more often than shifted versions of that same
pattern. Its p/q values are not posterior probabilities that a molecule exists.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

C_KMS = 299_792.458


def sigma_differences(values: np.ndarray) -> float:
    """Robust per-channel noise estimate from first differences (approximate)."""
    y = np.asarray(values, dtype=float)
    differences = np.diff(y[np.isfinite(y)])
    if differences.size < 3:
        return float("nan")
    center = np.median(differences)
    return float(1.4826 * np.median(np.abs(differences - center)) / np.sqrt(2))


def detect_features(
    frequency_mhz: np.ndarray,
    intensity: np.ndarray,
    *,
    threshold_sigma: float = 5.0,
    min_distance_channels: int = 3,
    absolute_prominence_floor: float = 0.03,
) -> tuple[pd.DataFrame, float, float]:
    """Find positive spectral peaks and return peaks, robust sigma, and channel step."""
    frequency = np.asarray(frequency_mhz, dtype=float)
    y = np.asarray(intensity, dtype=float)
    if frequency.size != y.size or frequency.size < 3:
        raise ValueError("Frecuencia e intensidad deben tener igual longitud y al menos 3 canales.")
    if not np.all(np.isfinite(frequency)) or not np.all(np.diff(frequency) > 0):
        raise ValueError("El eje de frecuencia debe ser finito y estrictamente creciente.")
    sigma = sigma_differences(y)
    step = float(np.median(np.diff(frequency)))
    if not np.isfinite(sigma) or sigma <= 0:
        return pd.DataFrame(columns=["canal", "frecuencia_pico_LSRK_MHz", "pico_Jy_beam",
                                     "prominencia_Jy_beam", "sigma_aprox_Jy_beam",
                                     "SNR_prominencia_aprox"]), sigma, step
    prominence = max(float(threshold_sigma) * sigma, float(absolute_prominence_floor))
    indices, properties = find_peaks(
        y, prominence=prominence, distance=max(1, int(min_distance_channels))
    )
    peaks = pd.DataFrame({
        "canal": indices,
        "frecuencia_pico_LSRK_MHz": frequency[indices],
        "pico_Jy_beam": y[indices],
        "prominencia_Jy_beam": properties["prominences"],
    })
    peaks["sigma_aprox_Jy_beam"] = sigma
    peaks["SNR_prominencia_aprox"] = peaks["prominencia_Jy_beam"] / sigma
    return peaks, sigma, step


def catalog_associations(
    catalog: pd.DataFrame,
    band_info: Mapping[int, tuple[float, float, float]],
    peaks_by_window: Mapping[int, pd.DataFrame],
    *,
    velocity_kms: float = 2.7,
    velocity_margin_kms: float = 1.0,
    min_snr: float = 5.0,
) -> pd.DataFrame:
    """Associate catalog transitions with nearest detected features.

    The frequency tolerance combines the stated velocity margin, laboratory
    uncertainty (when present), and half a channel. This is an operational
    association window, not a confidence interval.
    """
    rows: list[dict[str, object]] = []
    for window, (low, high, step) in band_info.items():
        peaks = peaks_by_window.get(window, pd.DataFrame())
        if peaks.empty:
            continue
        peak_freq = peaks["frecuencia_pico_LSRK_MHz"].to_numpy(dtype=float)
        for transition in catalog.itertuples(index=False):
            rest = float(transition.frecuencia_reposo_MHz)
            predicted = rest * (1 - float(velocity_kms) / C_KMS)
            if not low <= predicted <= high:
                continue
            index = int(np.argmin(np.abs(peak_freq - predicted)))
            peak = peaks.iloc[index]
            catalog_sigma = getattr(transition, "incertidumbre_MHz", np.nan)
            catalog_sigma = (float(catalog_sigma) if pd.notna(catalog_sigma)
                             and float(catalog_sigma) > 0 else 0.0)
            tolerance = np.sqrt(
                (predicted * float(velocity_margin_kms) / C_KMS) ** 2
                + catalog_sigma**2 + (step / 2) ** 2
            )
            tolerance = min(max(float(tolerance), step), 1.5)
            residual = float(peak_freq[index] - predicted)
            snr = float(peak["SNR_prominencia_aprox"])
            if abs(residual) <= tolerance and snr >= float(min_snr):
                channel = int(peak["canal"])
                rows.append({
                    "ventana": f"W{window}",
                    "especie": str(transition.especie),
                    "catalogo": str(transition.catalogo),
                    "entrada": str(transition.entrada),
                    "frecuencia_reposo_MHz": rest,
                    "frecuencia_predicha_LSRK_MHz": predicted,
                    "frecuencia_pico_LSRK_MHz": float(peak_freq[index]),
                    "residuo_MHz": residual,
                    "tolerancia_asociacion_MHz": tolerance,
                    "prominencia_Jy_beam": float(peak["prominencia_Jy_beam"]),
                    "SNR_aprox": snr,
                    "v_radio_kms": C_KMS * (rest - peak_freq[index]) / rest,
                    "grupo_rasgo": f"{window}-{channel}",
                })
    return pd.DataFrame(rows)


def _matched_feature_count(
    predicted: np.ndarray,
    tolerances: np.ndarray,
    peak_frequencies: np.ndarray,
    peak_channels: np.ndarray,
) -> int:
    if predicted.size == 0 or peak_frequencies.size == 0:
        return 0
    right = np.searchsorted(peak_frequencies, predicted)
    left_idx = np.clip(right - 1, 0, peak_frequencies.size - 1)
    right_idx = np.clip(right, 0, peak_frequencies.size - 1)
    nearest = np.where(
        np.abs(peak_frequencies[left_idx] - predicted)
        <= np.abs(peak_frequencies[right_idx] - predicted),
        left_idx, right_idx,
    )
    matched = np.abs(peak_frequencies[nearest] - predicted) <= tolerances
    return int(np.unique(peak_channels[nearest[matched]]).size)


def benjamini_yekutieli(p_values: Mapping[str, float]) -> dict[str, float]:
    """FDR-adjust p-values under arbitrary dependence (Benjamini–Yekutieli)."""
    if not p_values:
        return {}
    names = sorted(p_values, key=lambda key: p_values[key])
    count = len(names)
    harmonic = float(np.sum(1 / np.arange(1, count + 1)))
    adjusted: dict[str, float] = {}
    running = 1.0
    for rank in range(count, 0, -1):
        name = names[rank - 1]
        running = min(running, float(p_values[name]) * count * harmonic / rank, 1.0)
        adjusted[name] = running
    return adjusted


def assess_molecular_evidence(
    catalog: pd.DataFrame,
    band_info: Mapping[int, tuple[float, float, float]],
    peaks_by_window: Mapping[int, pd.DataFrame],
    *,
    velocity_kms: float = 2.7,
    velocity_margin_kms: float = 1.0,
    min_snr: float = 5.0,
    simulations: int = 2000,
    random_state: int = 20261011,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank catalog-alignment evidence using a circular-shift empirical null.

    The null shifts each species' transition pattern as a block within each
    covered window. The test controls for line density only approximately; it
    does not model line strengths, radiative transfer, or catalog completeness.
    Returned q-values are BY-adjusted empirical alignment p-values, never
    posterior probabilities of molecular identity.
    """
    species = sorted(catalog["especie"].dropna().astype(str).unique())
    associations = catalog_associations(
        catalog, band_info, peaks_by_window, velocity_kms=velocity_kms,
        velocity_margin_kms=velocity_margin_kms, min_snr=min_snr,
    )
    prepared: dict[str, dict[int, tuple[np.ndarray, np.ndarray]]] = {}
    for name in species:
        prepared[name] = {}
        for window, (low, high, step) in band_info.items():
            entries = catalog.loc[catalog["especie"].astype(str).eq(name)]
            rest = pd.to_numeric(entries["frecuencia_reposo_MHz"], errors="coerce").to_numpy()
            predicted = rest * (1 - float(velocity_kms) / C_KMS)
            in_band = (predicted >= low) & (predicted <= high) & np.isfinite(predicted)
            if not np.any(in_band):
                prepared[name][window] = (np.array([], dtype=float), np.array([], dtype=float))
                continue
            uncertainties = pd.to_numeric(entries["incertidumbre_MHz"], errors="coerce").to_numpy()
            uncertainties = np.where(np.isfinite(uncertainties) & (uncertainties > 0), uncertainties, 0.0)
            tolerance = np.sqrt(
                (predicted * float(velocity_margin_kms) / C_KMS) ** 2
                + uncertainties**2 + (step / 2) ** 2
            )
            tolerance = np.minimum(np.maximum(tolerance, step), 1.5)
            prepared[name][window] = (predicted[in_band], tolerance[in_band])

    observed: dict[str, int] = {}
    nulls = {name: np.zeros(max(1, int(simulations)), dtype=int) for name in species}
    peak_arrays: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for window in band_info:
        peaks = peaks_by_window.get(window, pd.DataFrame())
        peak_arrays[window] = (
            peaks.get("frecuencia_pico_LSRK_MHz", pd.Series(dtype=float)).to_numpy(dtype=float),
            peaks.get("canal", pd.Series(dtype=int)).to_numpy(dtype=int),
        )
    for name in species:
        observed[name] = sum(
            _matched_feature_count(*prepared[name][window], *peak_arrays[window])
            for window in band_info
        )

    rng = np.random.default_rng(random_state)
    n_sim = max(1, int(simulations))
    for replicate in range(n_sim):
        for name in species:
            total = 0
            for window, (low, high, _step) in band_info.items():
                predicted, tolerances = prepared[name][window]
                if predicted.size == 0:
                    continue
                offset = rng.uniform(0.0, high - low)
                shifted = low + np.mod(predicted - low + offset, high - low)
                total += _matched_feature_count(
                    shifted, tolerances, *peak_arrays[window]
                )
            nulls[name][replicate] = total

    empirical_p = {
        name: (1 + int(np.sum(nulls[name] >= observed[name]))) / (n_sim + 1)
        for name in species
    }
    adjusted_q = benjamini_yekutieli(empirical_p)
    groups_per_feature = (associations.groupby("grupo_rasgo")["especie"].nunique()
                          if not associations.empty else pd.Series(dtype=int))
    results: list[dict[str, object]] = []
    for name in species:
        current = associations.loc[associations["especie"].eq(name)].copy()
        if not current.empty:
            current["residuo_normalizado"] = (
                current["residuo_MHz"].abs() / current["tolerancia_asociacion_MHz"]
            )
            best = current.sort_values("residuo_normalizado").drop_duplicates("grupo_rasgo")
            velocities = best["v_radio_kms"].to_numpy(dtype=float)
            velocity_median = float(np.median(velocities))
            velocity_mad = float(1.4826 * np.median(np.abs(velocities - velocity_median)))
            shared_fraction = float(np.mean([
                groups_per_feature.get(group, 1) > 1 for group in best["grupo_rasgo"]
            ]))
            median_snr = float(best["SNR_aprox"].median())
        else:
            best = current
            velocity_median = np.nan
            velocity_mad = np.nan
            shared_fraction = np.nan
            median_snr = np.nan
        independent = int(observed[name])
        velocity_consistent = (
            np.isfinite(velocity_median)
            and abs(velocity_median - float(velocity_kms)) <= float(velocity_margin_kms)
            and np.isfinite(velocity_mad) and velocity_mad <= 0.5
        )
        if independent >= 3 and adjusted_q[name] <= 0.05 and velocity_consistent:
            grade = "Compatibilidad tentativa" if shared_fraction < 0.5 else "Sugestiva, ambigua"
        elif independent >= 2 and adjusted_q[name] <= 0.10 and velocity_consistent:
            grade = "Indicio débil; insuficiente para asignar"
        else:
            grade = "Sin evidencia suficiente en esta prueba"
        results.append({
            "especie": name,
            "rasgos_distintos_asociados": independent,
            "SNR_mediana_aprox": median_snr,
            "v_mediana_kms": velocity_median,
            "dispersion_v_MAD_kms": velocity_mad,
            "fraccion_rasgos_compartidos": shared_fraction,
            "nulo_mediana_rasgos": float(np.median(nulls[name])),
            "nulo_p95_rasgos": float(np.quantile(nulls[name], 0.95)),
            "p_empirico_desplazamiento": empirical_p[name],
            "q_BY_12_especies": adjusted_q[name],
            "grado_evidencia": grade,
            "interpretacion": "q ajustado evalúa alineación de frecuencias bajo el nulo circular; no es P(molécula|datos).",
        })
    summary = pd.DataFrame(results).sort_values(
        ["q_BY_12_especies", "rasgos_distintos_asociados"],
        ascending=[True, False], ignore_index=True,
    )
    return summary, associations


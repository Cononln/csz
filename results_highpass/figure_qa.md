# Q1 high-pass figures: render and data QA

- Delivered format: 11 PNG figures at 600 dpi. Figure-level PDF/SVG files were used only for temporary rendered checks and are not part of the delivery.
- Data coverage: 48 real-data summary rows (12 record–task–channel combinations × 4 cutoffs), 19,728 real-waveform rows (411 samples per combination/cutoff), 3,456 semi-synthetic recovery rows (12 replicates), and 6,168 example-waveform rows. All expected keys are present and unique. No source observations were dropped for display; the two A–Task 2–F3 plots and one simulated waveform are explicitly labeled examples.
- Numeric trace: the real-waveform Right − Left curve equals `right - left` at every analysis sample to numerical precision, and its 250–500 ms mean matches the summary table (maximum absolute discrepancy < 1e-8 source units). Simulation rows cover Left, Right, and Delta; plotted RMSE/correlation combine Left and Right, while difference-error plots use Delta. Each replicate is first summarized within its record–task–channel–condition combination, then across 24 response or 12 difference combinations; error bars are interquartile ranges across combinations. These combinations are descriptive strata, not independent biological replicates.
- Axes and units: amplitudes and errors remain in source units. Baseline RMS is explicitly labeled as residual −0.2 to 0 s RMS. Frequency response is the magnitude of the sequential two-pass high-pass and 30 Hz low-pass filters. The impulse plot shows absolute response on a log axis; values below 1e-7 are hidden only to avoid drawing a floor line.
- Source-code preflight: 18 PASS, 3 WARN, 0 FAIL. WARN review: TIFF is absent because PNG is the requested output; 172.7 mm is an analysis-figure width rather than a specified journal width; the logarithm in the frequency-response plot has an explicit positive floor (`np.maximum(magnitude, 1e-10)`).
- Rendered QA: all 11 temporary PDFs passed the 5 pt text audit (minimum observed 6.5 pt). All 11 passed the PDF geometry collision audit with zero FAIL and zero WARN. The panel-alignment gate records `NOT APPLICABLE` because each figure has one plot area; the heatmap colorbar is excluded from panel comparison. Compact `.alignment.json` and `.collision.json` checks remain beside each PNG.

| Figure | Visual inspection |
| --- | --- |
| `filter_frequency_response.png` | Four curves and −6 dB guide are legible; frequency and gain axes use the intended logarithmic/linear scales. |
| `filter_impulse_response.png` | Four absolute impulse tails are visible without a plotted floor or legend/tick collision. |
| `A_task2_F3_left_right.png` | Solid/dashed condition mapping, cue onset, and 250–500 ms interval are clear. |
| `A_task2_F3_delta.png` | Four Right − Left trajectories remain distinguishable across −0.2 to 0.8 s. |
| `delta_mean_full_matrix.png` | All 12 combinations and four cutoffs are present; signed values remain legible on the diverging scale. |
| `split_half_by_highpass.png` | Twelve paired trajectories and the median are legible without hiding negative correlations. |
| `baseline_noise_by_highpass.png` | Twelve paired trajectories and the median are legible; baseline interval is stated. |
| `simulation_recovery_example.png` | Reference, noisy pre-high-pass curve, and four recovered trajectories are distinct; this is a single example, not an aggregate. |
| `simulation_rmse.png` | Two reference families, center and IQR are visible; 24 response combinations underlie each summary. |
| `simulation_waveform_correlation.png` | Correlation contrast and IQR are visible without truncating the 0.5 Hz spread. |
| `simulation_delta_error.png` | Absolute condition-difference errors and IQR are visible; 12 combinations underlie each summary. |

These rendering checks do not establish biological ERP ground truth or choose a high-pass setting. The numerical interpretation belongs in the decision report.

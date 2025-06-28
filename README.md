# Virtual Cell Challenge
## Progres - 28 Juni 2025

### simple EDA
| Metrics | Nilai |
|--------|-------|
| Jumlah sel (training) | **221 273** |
| Jumlah gen (training) | **18 080** |
| Sel yang dipakai buat QC violin | **10 000** |

> Plot violin nunjukkin distribusi total UMI per sel (`n_counts`) & jumlah gen terdeteksi (`n_genes`) buat sampel random 10k sel.  
> Mayoritas sel ada di kisaran **25–40 k UMI counts** dan **8–9 k gen**, dengan bbrp sel berkualitas rendah di ujung distribusi.

## Rangkuman submission 1 (zero-shot, exploratory) 29 Juni 2025

### Strategi
Pake cell-mean doang
- grup traindata by `target_gene` (gen yang diperturb)
- hitung **mean expression profile** (buat 18,080 gen) utk tiap jenis perturbation
- Untuk tiap perturbation di validation, kita **assign mean profile** sbg prediksi awal

### Output
- Prepare submission pakai `cell-eval` dan pass the check

### File2 penting
- `baseline_mean.py`: script utk buat hitung baseline mean
- `pred_validation_meanexp_ntc.h5ad`: hasil prediksi dgn kontrol
- `my_first_submission.vcc`: hasil prediksi dgn format lomba
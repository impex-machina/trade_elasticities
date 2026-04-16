# Trade Elasticity Estimates

Replication code and documentation for heterogeneous trade elasticity estimates.
Data outputs are hosted on Hugging Face due to file size.

**Projects:**
- [trade_elast_baci_hs92_v202601_hs4](#trade_elast_baci_hs92_v202601_hs4) — Three-stage heterogeneous σ and γ estimator adapted for BACI at the country level, HS4 (Feenstra 1994; Soderbery 2018)
- [soderbery_2018_extension](#soderbery_2018_extension) — Direct extension of Soderbery (2018) joint nonlinear SUR estimator to BACI HS92 v202601, HS4

---

## `trade_elast_baci_hs92_v202601_hs4`

[![Dataset on HF](https://img.shields.io/badge/🤗%20Hugging%20Face-trade__elast__baci__hs92__v202601__hs4-yellow)](https://huggingface.co/datasets/impex-machina/trade_elast_baci_hs92_v202601_hs4)

Three-stage heterogeneous σ and γ estimator adapted for BACI at the country level, HS4, estimated at the country and regional level.

**Computation:** Scripts were run on AWS EC2 (`c7a.16xlarge`, us-east-1, 64 vCPUs, 128GB RAM). Full pipeline runtime is approximately 45 minutes on 62 cores. A 4-core machine is supported via serial fallback but expect 5–8 hours. Replicators will need to supply their own AMI (`[AMI-ID]`) and IAM role (`[IAM-ROLE]`). See the README in this folder for stage-level timing details and checkpoint/resume instructions.

| File | Description | Download |
|---|---|---|
| `baci_hs92_v202601_elast_country_hs4_fixed_sigma.rds` | Country-level elasticities, HS4, fixed σ | [↓](https://huggingface.co/datasets/impex-machina/trade_elast_baci_hs92_v202601_hs4/resolve/main/baci_hs92_v202601_elast_country_hs4_fixed_sigma.rds) |
| `baci_hs92_v202601_elast_regional_hs4_fixed_sigma.rds` | Regional elasticities, HS4, fixed σ | [↓](https://huggingface.co/datasets/impex-machina/trade_elast_baci_hs92_v202601_hs4/resolve/main/baci_hs92_v202601_elast_regional_hs4_fixed_sigma.rds) |

```r
tmp <- tempfile(fileext = ".rds")
download.file(
  "https://huggingface.co/datasets/impex-machina/trade_elast_baci_hs92_v202601_hs4/resolve/main/baci_hs92_v202601_elast_country_hs4_fixed_sigma.rds",
  destfile = tmp,
  mode = "wb"
)
df <- readRDS(tmp)
```

---

## `soderbery_2018_extension`

[![Dataset on HF](https://img.shields.io/badge/🤗%20Hugging%20Face-soderbery__2018__extension-yellow)](https://huggingface.co/datasets/impex-machina/soderbery_2018_extension)

Direct extension of the Soderbery (2018) joint nonlinear SUR estimator to BACI HS92 v202601, at the HS4 level, estimated at the country and regional level.

**Computation:** Scripts were run on AWS EC2 (`c7a.16xlarge`, us-east-1, 64 vCPUs, 128GB RAM). Replicators will need to supply their own AMI (`[AMI-ID]`) and IAM role (`[IAM-ROLE]`). See the README in this folder for runtime details and checkpoint/resume instructions.

| File | Description | Download |
|---|---|---|
| `baci_hs92_v202601_elast_country_hs4_step2.rds` | Country-level elasticities, HS4, step 2 | [↓](https://huggingface.co/datasets/impex-machina/soderbery_2018_extension/resolve/main/baci_hs92_v202601_elast_country_hs4_step2.rds) |
| `baci_hs92_v202601_elast_regional_hs4_step2.rds` | Regional elasticities, HS4, step 2 | [↓](https://huggingface.co/datasets/impex-machina/soderbery_2018_extension/resolve/main/baci_hs92_v202601_elast_regional_hs4_step2.rds) |

```r
tmp <- tempfile(fileext = ".rds")
download.file(
  "https://huggingface.co/datasets/impex-machina/soderbery_2018_extension/resolve/main/baci_hs92_v202601_elast_country_hs4_step2.rds",
  destfile = tmp,
  mode = "wb"
)
df <- readRDS(tmp)
```

---

## License

Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — free to use with attribution.  
Code: MIT

## Citation

If you use these estimates, please cite:
- Feenstra, R.C. (1994). New product varieties and the measurement of international prices. *American Economic Review*, 84(1), 157–177.
- Soderbery, A. (2018). Trade elasticities, heterogeneity, and optimal tariffs. *Journal of International Economics*, 114, 44–62.
- CEPII BACI World Trade Database, HS92 v202601.

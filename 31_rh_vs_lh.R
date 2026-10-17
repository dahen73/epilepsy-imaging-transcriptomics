# =====================================================================
# 31_rh_vs_lh.R -- 节点2nd: lh vs rh 对比汇总 (2026-10-15)
# 依赖: 30 号 rh 结果 + lh 已有结果 (10/16/19 号)
# 输出: full_run/rh/rh_vs_lh_summary.txt + rh_vs_lh.csv
# =====================================================================
suppressPackageStartupMessages({ library(data.table) })
ws     <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
rdir   <- file.path(ws, "full_run/rh")
lh_map <- list(                       # 模态: lh 结果文件前缀
  thick = list(dir = "full_run/targeted",      files = c("individual_targeted")),
  area  = list(dir = "full_run/targeted_area", files = c("individual_targeted")),
  vol   = list(dir = "full_run/targeted_vol",  files = c("individual_targeted")),
  sulc  = list(dir = "full_run/targeted_sulc", files = c("individual_targeted"))
)
lh_files <- list(
  thick = c(Number_ASMs = "Number_ASMs", age_onset = "age_onset", HS_vs_other = "HS_vs_other",
            SE = "SE", ILAE_Y1 = "ILAE_Y1", FCDII_vs_other = "FCDII_vs_other"),
  area  = c(Number_ASMs = "Number_ASMs", age_onset = "age_onset", HS_vs_other = "HS_vs_other",
            SE = "SE", ILAE_Y1 = "ILAE_Y1", FCDII_vs_other = "FCDII_vs_other"),
  vol   = c(Number_ASMs = "Number_ASMs", age_onset = "age_onset", HS_vs_other = "HS_vs_other",
            SE = "SE", ILAE_Y1 = "ILAE_Y1", FCDII_vs_other = "FCDII_vs_other"),
  sulc  = c(Number_ASMs = "Number_ASMs", age_onset = "age_onset", HS_vs_other = "HS_vs_other",
            SE = "SE", ILAE_Y1 = "ILAE_Y1", FCDII_vs_other = "FCDII_vs_other")
)
rows <- list()
for (mn in names(lh_map)) {
  for (pn in names(lh_files[[mn]])) {
    frh <- file.path(rdir, sprintf("individual_rh_%s_%s.csv", mn, pn))
    flh <- file.path(ws, lh_map[[mn]]$dir, sprintf("individual_targeted_%s.csv", pn))
    if (!file.exists(frh) || !file.exists(flh)) next
    d_rh <- fread(frh); d_lh <- fread(flh)
    m <- merge(d_rh[, .(gene_set, description, beta_rh = beta, p_emp_rh = p_emp, fdr_rh = p_emp_fdr)],
               d_lh[, .(gene_set, beta_lh = beta, p_emp_lh = p_emp, fdr_lh = p_emp_fdr)],
               by = "gene_set")
    m[, `:=`(mod = mn, pred = pn,
             dir_agree = sign(beta_rh) == sign(beta_lh),
             rh_fdr05 = fdr_rh < 0.05, lh_fdr05 = fdr_lh < 0.05)]
    rows[[length(rows) + 1]] <- m
  }
}
all <- rbindlist(rows)
fwrite(all, file.path(rdir, "rh_vs_lh.csv"))
lines <- c("================================================================",
           "右半球敏感性分析: lh (34 区) vs rh (32 区) 对比汇总",
           "================================================================", "")
# 全景矩阵
lines <- c(lines, "== 1. 全景: 各 mod x pred 的 FDR<0.05 集数 (lh vs rh) ==")
piv <- all[, .(lh_fdr = sum(lh_fdr05), rh_fdr = sum(rh_fdr05),
               lh_p05 = sum(p_emp_lh < 0.05), rh_p05 = sum(p_emp_rh < 0.05),
               dir_agree = sum(dir_agree), n = .N), by = .(mod, pred)]
print(piv)
lines <- c(lines, capture.output(print(as.data.frame(piv))), "")
# lh 显著集在 rh 中的表现
lines <- c(lines, "== 2. lh FDR 显著集在 rh 中的表现 ==")
sig <- all[lh_fdr05 == TRUE]
if (nrow(sig)) {
  agg <- sig[, .(dir_agree = sum(dir_agree), n = .N,
                 rh_fdr_sig = sum(rh_fdr05), rh_nominal = sum(p_emp_rh < 0.05),
                 med_r_fdr = round(median(fdr_rh), 3)), by = .(mod, pred)]
  print(as.data.frame(agg))
  lines <- c(lines, capture.output(print(as.data.frame(agg))), "")
  lines <- c(lines, sprintf("合计: lh 显著集 %d 个 | rh 方向一致 %d (%.0f%%) | rh FDR 显著 %d | rh 名义显著 %d",
                            nrow(sig), sum(sig$dir_agree), 100 * mean(sig$dir_agree),
                            sum(sig$rh_fdr05), sum(sig$p_emp_rh < 0.05)), "")
}
lines <- c(lines, "== 3. rh 新增 FDR 显著集 (lh 阴性) ==")
new_sig <- all[!lh_fdr05 & rh_fdr05]
if (nrow(new_sig)) { print(as.data.frame(new_sig[, .(mod, pred, gene_set, description, p_emp_rh, fdr_rh)])) } else { lines <- c(lines, "(无)") }
lines <- c(lines, "", "================================================================")
writeLines(lines, file.path(rdir, "rh_vs_lh_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")

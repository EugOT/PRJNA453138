library(tidyverse)
metrics_summary_SRR7059018 <- read_csv("/data/PRJNA453138/cellranger/SRR7059018/outs/metrics_summary.csv") %>% mutate(Run = "SRR7059018")
metrics_summary_SRR7059019 <- read_csv("/data/PRJNA453138/cellranger/SRR7059019/outs/metrics_summary.csv") %>% mutate(Run = "SRR7059019")
metrics_summary_SRR7059020 <- read_csv("/data/PRJNA453138/cellranger/SRR7059020/outs/metrics_summary.csv") %>% mutate(Run = "SRR7059020")
metrics_summary_SRR7059021 <- read_csv("/data/PRJNA453138/cellranger/SRR7059021/outs/metrics_summary.csv") %>% mutate(Run = "SRR7059021")
metrics_summary_SRR7059022 <- read_csv("/data/PRJNA453138/cellranger/SRR7059022/outs/metrics_summary.csv") %>% mutate(Run = "SRR7059022")
metrics_summary_SRR7059023 <- read_csv("/data/PRJNA453138/cellranger/SRR7059023/outs/metrics_summary.csv") %>% mutate(Run = "SRR7059023")
metrics_summary <-
    bind_rows(
        metrics_summary_SRR7059018,
        metrics_summary_SRR7059019,
        metrics_summary_SRR7059020,
        metrics_summary_SRR7059021,
        metrics_summary_SRR7059022,
        metrics_summary_SRR7059023)

metrics_summary |>
    select("Estimated Number of Cells", "Run")

write_tsv(metrics_summary, here::here("metrics_summary.tsv"))


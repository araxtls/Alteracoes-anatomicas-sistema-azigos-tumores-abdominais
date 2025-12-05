# CÁLCULO DE TAMANHO MÍNIMO DE AMOSTRA
valor_critico = 1.96 ##para z score de 0.95
desvio1_mm = 1
desvio2_mm = 2
erro_mm = 0.5

N_criancas_menores = (valor_critico*desvio1_mm / erro_mm)^2 %>%
  ceiling()
N_criancas_menores

N_criancas_maiores = (valor_critico*desvio2_mm / erro_mm)^2 %>%
  ceiling()
N_criancas_maiores

# ESTATÍSTICA DESCRITIVA
library(tidyverse)
dados = read.csv("dados_pesquisa1.csv")
dados <- dados %>%
  mutate(
    faixa_etaria = cut(
      idade,
      breaks = seq(0, 12, by = 2),
      include.lowest = TRUE,
      right = FALSE,
      labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12")
    )
  ) ## estratificando os pacientes em grupos de idade

cont_vars <- c("idade_anos","maior_diam_tumor_mm","diam_VCS_mm","den_VCS_hu","dist_vert_VCS_mm",
               "diam_VCI_mm","den_VCI_hu","dist_vert_VCI_mm","diam_azigo_mm","den_azigo_hu",
               "dist_vert_azigo_mm","hemiazigo_diam_mm","den_hemiazigos_hu","dist_vert_hemiazigo_mm",
               "diam_tributaria_mm","qtde_tributarias","diam_vcs_mm","den_vcs_hu","diam_vci_mm","den_vci_hu") ##selecionando variaveis continuas
cat_vars <- c("sexo","tipo_tumor","localizacao_tumor","vis_tributarias","vis_hemiazigos") #selecionando variaveis categoricas

#1 Descricao variaveis categoricas
freq_table <- dados %>%
  select(all_of(cat_vars)) %>%
  pivot_longer(every_of(cat_vars), names_to = "variable", values_to = "value") %>%
  group_by(variable, value) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n)) ##tabela de frequencia
write_csv(freq_table, file.path(dir_out, "freq_table.csv"))

cat_comparisons <- list()
for(v in cat_vars){
  tab <- table(df$faixa_idade, df[[v]]) ##comparando variaveis categoricas DE ACORDO COM ESTRATOS DE IDADE
  test <- tryCatch({
    if(any(tab < 5)) ##verificando normalidade da distribuicao
      fisher.test(tab)
    else chisq.test(tab)
  }, error = function(e) list(p.value = NA))
  cat_comparisons[[v]] <- tibble(variable = v, p_value = test$p.value)
}
cat_comp_df <- bind_rows(cat_comparisons)
write_csv(cat_comp_df, file.path(dir_out, "comparisons_categorical.csv"))

#2 Descricao variaveis continuas
desc_overall <- dados %>%
  select(any_of(cont_vars)) %>%
  pivot_longer(every_of(intersect(cont_vars, names(df))),
               names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    n = sum(!is.na(value)),
    media = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) ##medidas de dispersao e centralidade GLOBAIS
write_csv(desc_overall, file.path(dir_out, "desc_overall.csv"))

desc_by_strata <- dados %>%
  select(faixa_idade, any_of(cont_vars)) %>%
  pivot_longer(cols = -faixa_idade, names_to = "variable", values_to = "value") %>%
  group_by(faixa_idade, variable) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) ##medidas de dispersao e centralidade DE ACORDO COM OS ESTRATOS DE IDADE
write_csv(desc_by_strata, file.path(dir_out, "desc_by_strata.csv"))

comparisons <- list()
for(var in intersect(cont_vars, names(dados))){ ##comparando variaveis categoricas DE ACORDO COM ESTRATOS DE IDADE
  sub <- dados %>%
    select(faixa_idade, !!sym(var)) %>%
    filter(!is.na(!!sym(var)))
  counts <- sub %>%
    group_by(faixa_idade) %>%
    tally() %>% pull(n)
  if(length(counts) < 2 || any(counts < 3)){
    comparisons[[var]] <- tibble(variable = var, method = NA_character_, p_value = NA_real_)
    next
  }

    shapiro_p <- sub %>% ##testando normalidade
      group_by(faixa_idade) %>%
    summarise(p = shapiro_test(!!sym(var))$p.value, .groups = "drop")
  use_t_test <- all(shapiro_p$p > 0.05, na.rm = TRUE)
  if(use_t_test){
    res <- t_test(sub, formula = as.formula(paste(var, "~ faixa_idade")), var.equal = FALSE)
    comparisons[[var]] <- res %>% select(statistic, p) %>% mutate(variable = var, method = "t_test")
  } else {
    res <- wilcox_test(sub, formula = as.formula(paste(var, "~ faixa_idade")))
    comparisons[[var]] <- res %>% select(statistic, p) %>% mutate(variable = var, method = "wilcox")
  }
}
comparisons_df <- bind_rows(comparisons)
write_csv(comparisons_df, file.path(dir_out, "comparisons_continuous.csv"))
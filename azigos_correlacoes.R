install.packages("rstatix")
library(tidyverse)
library(rstatix)

df <- read.csv("dados_pesquisa1.csv", stringsAsFactors = FALSE)
diametros <- c("diam_azigos_mm","hemiazigo_diam_mm","diam_cava_superior_mm",
               "diam_cava_inferior_mm","diam_tributaria_mm")

densidades <- c("den_azigos_hu","den_hemiazigos_hu","den_cava_superior_hu",
                "den_cava_inferior_hu","den_tributaria_hu")

outras_vars <- c("maior_diam_tumor_mm", "idade_anos")

calc_cor <- function(x, y, data){
  tryCatch({
    res <- cor.test(data[[x]], data[[y]], method = "spearman", exact = FALSE)
    tibble(var_x = x, var_y = y,
           rho = res$estimate,
           p_value = res$p.value)
  }, error = function(e){
    tibble(var_x = x, var_y = y, rho = NA, p_value = NA)
  })
}
calc_cor() ##correlacao de spearman

# 1 Correlacionando diametro dos vasos
cor_vasos_diam <- expand.grid(diam_x = diametros, diam_y = diametros) %>%
  filter(diam_x != diam_y) %>%
  pmap_dfr(~calc_cor(..1, ..2, df))
cor_vasos_diam

# 2 Correlacionando densidade dos vasos
cor_vasos_den <- expand.grid(den_x = densidades, den_y = densidades) %>%
  filter(den_x != den_y) %>%
  pmap_dfr(~calc_cor(..1, ..2, df))
cor_vasos_den

# 3 Correlacionando diametro e densidade no mesmo vaso
pares_vasos <- tibble(
  diam = diametros,
  dens = densidades
)

cor_diam_den <- pmap_dfr(pares_vasos, ~calc_cor(..1, ..2, df))
cor_diam_den

# 4 Diametro dos vasos x tamanho do tumor
cor_diam_tumor <- map_dfr(diametros, ~calc_cor(.x, "maior_diam_tumor_mm", df))
cor_diam_tumor

# 5 Densidade dos vasos x tamanho do tumor
cor_den_tumor <- map_dfr(densidades, ~calc_cor(.x, "maior_diam_tumor_mm", df))
cor_den_tumor

# 6 Diametro dos vasos x idade do paciente
cor_diam_idade <- map_dfr(diametros, ~calc_cor(.x, "idade_anos", df))
cor_diam_idade

# 7 Densidade dos vasos x tamanho do tumor
cor_den_idade <- map_dfr(densidades, ~calc_cor(.x, "idade_anos", df))
cor_den_idade
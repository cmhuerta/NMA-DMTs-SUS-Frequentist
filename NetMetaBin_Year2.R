############################################################################

# Metanálise em Rede Frequentista de TMDs disponíveis no SUS
#
# Por:   Caio Huerta
# Data: 2024-01-16
# última atualização: 2024-01-29


############################################################################


# instala pacotes necessários para a análise

install.packages("netmeta")
install.packages("rgl")
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("MathiasHarrer/dmetar")
install.packages("writexl")


############################################################################


# carrega pacotes 

library(netmeta)
library(rgl)
library(dmetar)
library(writexl)


############################################################################

# importa base de dados em formato .CSV

meta_year2 <- read.csv("Analysis_Year2.csv", sep = ";")


############################################################################

# prepara base para condução de metanálise


# adiciona legenda dos tratamentos
long.labels <- c("alemtuzumab", "fingolimode", "fumarato", "glatiramer", "IFNB1a_22", "IFNB1a_30", "IFNB1a_44", "IFNB1b",
                 "natalizumabe", "placebo", "teriflunomida14", "teriflunomida7")


# transforma dados do formato arm-based para formato de contraste
# sm = "OR" define odds ratio como medida sumário
# leganda: TE =  log odds ratio, seTE = standard errors  
p2 <- pairwise(treat = treatment, event = event.free, n = randomized, studlab = study,
               data = meta_year2, sm = "OR")


# netmetabin() solicita metanalise de desfecho binário
# calcula modelo de efeito fixo e aleatório
# método de inverso da variância (method = "Inverse)
# placebo como referência para as comparações pareadas (ref = "placebo")
nb2 <- netmetabin(p2, ref = "placebo", method = "Inverse")

############################################################################

# resultados


# base de de dados no formato de contraste
# coluna TE apresenta efeito de tratamento em log odds ratio
# coluna seTE apresenta erro padrão
p2


# sumário das análises de efeito fixo e efeitos aleatórios
nb2
summary(nb2)


# league table 
# triângulo inferior apresenta estimativa de evidência direta e indireta
# triângulo superior apresenta estimativa de evidência direta 
netleague <- netleague(nb2, 
                       bracket = "(", # use round brackets
                       digits=2)      # round to two digits
write.csv(netleague$fixed, "csv/netleaguefixed_ano2.csv")
write.csv(netleague$random, "csv/netleaguerandom_ano2.csv")

# código para salvar direto em xlsx com o pacote "writexl"
netleague(nb2, path = "netleaguetest2.xlsx")


# rede de evidência
jpeg(filename = "images/network_ano2.jpg", width = 800, height = 600, quality = 100)
netgraph(nb2, 
         number = TRUE,
         labels = paste0(long.labels, "\n(n=", round(n.trts), ")"),
         points = TRUE,
         cex.points = n.trts,
         offset = 0.05)
invisible(dev.off())


# rede de evidência em 3 dimensões
netgraph(nb2, dim = "3d")


# forest plot, common effect
jpeg(filename = "images/forest_common_ano2.jpg", width = 600, height = 400, quality = 100)
forest(nb2,
       reference.group = "placebo",
       sortvar = TE,
       xlim = c(0.1, 10),
       smlab = paste("TMDs vs. placebo \n",
                     "Common effect model"),
       drop.reference.group = TRUE,
       leftcols = c("studlab","k"), 
       leftlabs = c("Drug","Direct\nstudies"),
       label.left = "Favors placebo",
       label.right = "Favors Intervention",
       labels = long.labels)
invisible(dev.off())

# forest plot, random-effects 
jpeg(filename = "images/forest_random_ano2.jpg", width = 600, height = 400, quality = 100)
forest(nb2,
       reference.group = "placebo",
       pooled = "random",
       sortvar = TE,
       xlim = c(0.1, 10),
       smlab = paste("TMDs vs. placebo \n",
                     "Random-effects model"),
       drop.reference.group = TRUE,
       leftcols = c("studlab","k"), 
       leftlabs = c("Drug","Direct\nstudies"),
       label.left = "Favors placebo",
       label.right = "Favors Intervention",
       labels = long.labels)
invisible(dev.off())

# ranqueamento, método SUCRA
netrankSUCRA <- netrank(nb2, method = "SUCRA", small.values = "bad",
                        bracket = "(")
write.csv(netrankSUCRA$ranking.fixed, "csv/SUCRAnetrankfixed_ano2.csv")
write.csv(netrankSUCRA$ranking.random, "csv/SUCRAnetrankrandom_ano2.csv")

# ranqueamento, método P-score
netrankP <- netrank(nb2, method = "P-score", small.values = "bad",
                    bracket = "(")
write.csv(netrankP$ranking.fixed, "csv/Pnetrankfixed_ano2.csv")
write.csv(netrankP$ranking.random, "csv/Pnetrankrandom_ano2.csv")

# rankograma

ran2 <- rankogram(nb2, nsim = 1000,small.values = "undesirable")
print(ran2, cumulative.rankprob = TRUE)
plot(ran2)

# gráfico netspliting para avaliar inconsistência da rede
# estima a contribuição da evidência direta e indireta de cada contraste com evidência direta disponível
jpeg("images/netsplitting_ano2.jpg", width = 900, height = 1200, quality = 100)
forest(netsplit(nb2))
invisible(dev.off())

# gráfico avalia proporção de evidência direta de todos os contrastes possíveis
jpeg("images/direct_ano2.jpg", width = 900, height = 600, quality = 100)
d.evidence2 <- direct.evidence.plot(nb2)
plot(d.evidence2)
invisible(dev.off())

#decomp.design
decomp.design(nb2)


############################################################################
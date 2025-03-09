############################################################################

# Metanálise em Rede Frequentista de TMDs disponíveis no SUS
#
# Por:   Caio Huerta


############################################################################

# instala pacotes necessários para a análise

install.packages("netmeta")
install.packages("rgl")
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("MathiasHarrer/dmetar")

############################################################################

# carrega pacotes 

library(netmeta)
library(rgl)
library(dmetar)

############################################################################

# importa base de dados em formato .CSV
# base estruturada com as seguintes colunas para NMA de dados binários: study, treatment, event-free, randomized
db.frequentist <- read.csv("Database_Frequentist.csv", sep = ";")

############################################################################

# preparação da base de dados para condução de metanálise


# adiciona legenda dos tratamentos
labels <- c("GLA", "ALE", "FIN", "FUM", "IFN-1a IM", "IFN-1a SC 22", "IFN-1a SC 44",
            "IFN-1b", "NAT", "PBO", "TER14", "TER7")
 
# transforma base de dados do formato arm-based (1 linha contém informação de cada braço do estudo) para formato de contraste (1 linha para cada comparação possível do estudo)
pairwise.freq <- pairwise(treat = treatment, 
                   event = event.free, 
                   n = randomized, studlab = study,
                   data = db.frequentist, sm = "OR")  # sm = "OR" define odds ratio como medida sumário 

# netmetabin() solicita metanalise de desfecho binário
# calcula modelo de efeito fixo e aleatório
# método de inverso da variância (method = "Inverse)
# placebo como referência para as comparações pareadas (ref = "placebo")
nb <- netmetabin(pairwise.freq, ref = "placebo", method = "Inverse")
nb.common <- netmetabin(pairwise.freq, ref = "placebo", method = "Inverse", random = FALSE)

############################################################################

# resultados

# visualiza base de dados no formato de contraste
# leganda: TE =  log odds ratio, seTE = standard errors 
View(pairwise.freq)


# sumário das análises de efeito fixo e efeitos aleatórios
summary(nb)

# rede de evidência
jpeg(filename = "images/network.jpg", width = 800, height = 600, quality = 100)
netgraph(nb, 
         number = TRUE,
         labels = paste0(labels, "\n(n=", round(n.trts), ")"),
         points = TRUE,
         cex.points = n.trts,
         col.points = "mediumblue", 
         lwd = 3, 
         plastic = FALSE,
         offset = 0.04) 
invisible(dev.off())

# rede de evidência em 3 dimensões
netgraph(nb, dim = "3d")


# league table 
# triângulo inferior apresenta estimativa de evidência direta e indireta
# triângulo superior apresenta estimativa de evidência direta 
options(digits = 2, OutDec = ",")
netleague <- netleague(nb, 
                       bracket = "(", # use round brackets
                       digits=2)      # round to two digits
write.csv(netleague$fixed, "csv/netleaguefixed.csv")
write.csv(netleague$random, "csv/netleaguerandom.csv")

# código para salvar direto em xlsx com o pacote "writexl"
options(digits = 2, OutDec = ",")
netleague(nb, bracket = "(", digits = 2, path = "csv/netleaguetest2.xlsx", overwrite = TRUE)

# forest plot, common effect
jpeg(filename = "images/forest_common.jpg", width = 600, height = 400, quality = 100)
forest(nb,
       reference.group = "placebo",
       pooled = "common",
       sortvar = TE,
       xlim = c(0.1, 10),
       smlab = paste("TMDs vs. placebo \n",
                     "Common effect model"),
       drop.reference.group = TRUE,
       leftcols = c("studlab","k"), 
       leftlabs = c("Drug","Direct\nstudies"),
       label.left = "Favors placebo",
       label.right = "Favors Intervention")
invisible(dev.off())

# forest plot, random-effects 
jpeg(filename = "images/forest_random.jpg", width = 600, height = 400, quality = 100)
forest(nb,
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
       label.right = "Favors Intervention")
invisible(dev.off())

# ranqueamento, método SUCRA
options(OutDec = ",")
netrankSUCRA <- netrank(nb, method = "SUCRA", small.values = "bad",
                        bracket = "(")
write.csv(netrankSUCRA$ranking.fixed, "csv/SUCRAnetrankfixed.csv")
write.csv(netrankSUCRA$ranking.random, "csv/SUCRAnetrankrandom.csv")


# ranqueamento, método P-score
netrankP <- netrank(nb, method = "P-score", small.values = "bad",
                    bracket = "(")

# rankograma

jpeg(filename = "images/rankogram.fixed.jpg", width = 600, height = 400, quality = 100)
ran <- rankogram(nb, nsim = 1000,small.values = "undesirable", cumulative = FALSE, common = TRUE)
print(ran, cumulative.rankprob = FALSE)
plot(ran)
invisible(dev.off())

jpeg(filename = "images/rankogram.random.jpg", width = 600, height = 400, quality = 100)
ran <- rankogram(nb, nsim = 1000,small.values = "undesirable", cumulative = FALSE, random = TRUE)
print(ran, cumulative.rankprob = FALSE)
plot(ran)
invisible(dev.off())

# gráfico netspliting para avaliar inconsistência da rede
# estima a contribuição da evidência direta e indireta de cada contraste com evidência direta disponível
jpeg("images/netsplitting_random.jpg", width = 900, height = 1200, quality = 100)
forest(netsplit(nb))
invisible(dev.off())

jpeg("images/netsplitting_common.jpg", width = 900, height = 1200, quality = 100)
forest(netsplit(nb.common))
invisible(dev.off())


netsplit(nb)
netsplit(nb.common)

# gráfico avalia proporção de evidência direta de todos os contrastes possíveis

jpeg("images/direct.evidence_common.jpg", width = 1800, height = 1200, quality = 90)
d.evidence.common <- direct.evidence.plot(nb)
plot(d.evidence.common)
invisible(dev.off())

jpeg("images/direct.evidence_random.jpg", width = 1800, height = 1200, quality = 90)
d.evidence.random <- direct.evidence.plot(nb, random = TRUE)
plot(d.evidence.random)
invisible(dev.off())


############################################################################
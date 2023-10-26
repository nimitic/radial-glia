library(reshape2)
library(ggplot2)

## re-making 4B pie charts as bar plots.
# copying the numbers from Nina's python notebook (scv_tel_clust_figs.ipynb) to avoid setting up the kernel...



transition_probabilities <- data.frame(c(0.9846314902992859, 0.011828193733383565, 8.922597778008176e-05, 0.0, 0.0, 0.00033523619942823574, 0.003115853790122243),
           c(0.24021123855843743, 0.7547824520184409, 0.003653247930704281, 0.0, 0.0013530614924174335, 0.0, 0.0),
           c(0.047170900111847304, 0.08974324813836333, 0.5827913628000015, 0.0, 0.10884206383310861, 0.14659845313564493, 0.024853971981034235),
           c(0.017004648565963856, 0.010970066468570717, 0.019449980508197906, 0.5908398017959375, 0.010838167134460164, 0.24484609125105855, 0.10605124427581132),
           c(0.03811058307563654, 0.11361745527475438, 0.08807324181377414, 0.0040088660115795835, 0.44499357546245744, 0.29163374652915486, 0.01956253183264309),
           c(0.013682263908182253, 0.0074388563086144315, 0.004835339923793374, 0.07224188008940696, 0.028116090901555018, 0.8531993296272505, 0.020486239241197374),
           c(0.04936373101730975, 0.009920685968780292, 0.021487366815760123, 0.060267904027412396, 0.007776922916275916, 0.07601839691974621, 0.7751649923347153),
           row.names = c('Neurons (other)', 'Neurons newborn', 'Proliferating cells',
                         'Radial glia (other)', 'Radial glia gfap++', 'Radial glia id2b+',
                         'Radial glia snap25a+'))

colnames(transition_probabilities) <- c('Neurons (other)', 'Neurons newborn', 'Proliferating cells',
                                        'Radial glia (other)', 'Radial glia gfap++', 'Radial glia id2b+',
                                        'Radial glia snap25a+')

transition_probabilities$TransitionTo <- rownames(transition_probabilities)

transition_probabilities <- melt(transition_probabilities, ID = "TransitionTo")

# rename gafp++ to her4++
transition_probabilities$variable <- gsub(transition_probabilities$variable, pattern = "Radial glia gfap++", replacement = "Radial glia her4++", fixed = TRUE)
transition_probabilities$TransitionTo <- gsub(transition_probabilities$TransitionTo, pattern = "Radial glia gfap++", replacement = "Radial glia her4++", fixed = TRUE)

## order like in curent figure
transition_probabilities$variable <- factor(transition_probabilities$variable, levels = c('Neurons (other)', 'Neurons newborn', 'Proliferating cells',
                                                       'Radial glia her4++', 'Radial glia id2b+', 'Radial glia (other)',
                                                       'Radial glia snap25a+'))
transition_probabilities$TransitionTo <- factor(transition_probabilities$TransitionTo, levels = c('Neurons (other)', 'Neurons newborn', 'Proliferating cells',
                                                                                          'Radial glia her4++', 'Radial glia id2b+', 'Radial glia (other)',
                                                                                          'Radial glia snap25a+'))

celltype_colours <- c('Neurons (other)' = "#ffffcc",
                      'Neurons newborn' = "#ffff33",
                      'Proliferating cells' = "#ff7f00",
                      'Radial glia her4++' = "#984ea3",
                      'Radial glia id2b+' = "#377eb8",
                      'Radial glia (other)' = "#bababa",
                      'Radial glia snap25a+' = "#e41a1c")

ggplot(data = transition_probabilities, mapping = aes(x = variable, y = value, fill = TransitionTo)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = celltype_colours, name = "transition to") +
  theme_minimal() +
  xlab("transition from") +
  ylab("fraction of cells")

pdf(file = "20230727_RG_4B_barcharts.pdf", width = 11, height = 3)
last_plot()
dev.off()

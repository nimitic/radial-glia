library(reshape2)
library(ggplot2)

## re-making 4B pie charts as bar plots.
# copying the numbers from Nina's python notebook (scv_tel_clust_figs.ipynb) to easier match plotting style in R



transition_probabilities <- data.frame(
           c(0.9857994786825948, 0.009741895478346694, 0.0004919811436818811, 0.0, 0.0, 0.0, 0.003966644695376653),
           c(0.22471965031666857, 0.77356308592672, 0.0, 0.0, 0.0, 0.0, 0.0017172637566114737),
           c(0.05826780693845225, 0.09556284393583657, 0.6092259407616354, 0.005854126353607041, 0.07244942746482981, 0.12799255788070332, 0.030647296664935635),
           c(0.014715707179188652, 0.011415924658196686, 0.0131685643823175, 0.5842661949845215, 0.00539226415180364, 0.24590977397200683, 0.12513157067196518),
           c(0.03971977130176595, 0.06628844535233099, 0.12909179867552475, 0.0, 0.4660559389460669, 0.2712398551407272, 0.02760419058358425),
           c(0.007724742774431791, 0.006517702979240734, 0.005680302438384393, 0.05712635140605482, 0.032408348532739294, 0.8678152907454717, 0.022727261123677275),
           c(0.04821255454000896, 0.013960751595219189, 0.016668546559538866, 0.06990914602841211, 0.014281436640954766, 0.0762183265258197, 0.7607492381100465),
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

## order like in current figure
transition_probabilities$variable <- factor(transition_probabilities$variable, levels = c('Neurons (other)', 'Neurons newborn', 'Proliferating cells',
                                                       'Radial glia her4++', 'Radial glia id2b+', 'Radial glia (other)',
                                                       'Radial glia snap25a+'))
transition_probabilities$TransitionTo <- factor(transition_probabilities$TransitionTo, levels = c('Neurons (other)', 'Neurons newborn', 'Proliferating cells',
                                                                                          'Radial glia her4++', 'Radial glia id2b+', 'Radial glia (other)',
                                                                                          'Radial glia snap25a+'))

celltype_colours <- c('Neurons (other)' = "#e69f00",
                      'Neurons newborn' = "#d55e00",
                      'Proliferating cells' = "#000000",
                      'Radial glia her4++' = "#cc79a7",
                      'Radial glia id2b+' = "#009e73",
                      'Radial glia (other)' = "#bbbbbb",
                      'Radial glia snap25a+' = "#0072b2")

ggplot(data = transition_probabilities, mapping = aes(x = variable, y = value, fill = TransitionTo)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = celltype_colours, name = "transition to") +
  theme_minimal() +
  xlab("transition from") +
  ylab("fraction of cells") +
  theme(text = element_text(size=14), legend.text=element_text(size=14), axis.text.x = element_text(angle = 10, vjust = 1, hjust=1))

pdf(file = "20230727_RG_4B_barcharts.pdf", width = 11, height = 3)
last_plot()
dev.off()

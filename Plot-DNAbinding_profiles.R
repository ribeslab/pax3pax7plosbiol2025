install.packages("ggseqlogo")

library(ggseqlogo)
library(ggplot2)

#A C G T
data(ggseqlogo_sample)
pfm <- matrix(c(
  0.01	, 0.01	, 0.97	, 0.01
  , 0.01	, 0.01	, 0.97	, 0.01
  , 0.01	, 0.97	, 0.01	, 0.01
  , 0.01	, 0.48	, 0.48	, 0.01
  , 0.01	, 0.48	, 0.48	, 0.01
  , 0.25	, 0.25	, 0.25	, 0.25
), ncol = 4, byrow = TRUE)

pfm <- t(pfm)
rownames(pfm) <- c("A", "C", "G", "T")

seqlogo_obj <- ggseqlogo(pfm, method = "bits")

seqlogo_obj



pdf(file = "Smad1-5-8-Martin-Malpartida2017.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 5)

seqlogo_obj


dev.off()

------------------------------
seqlogo_obj + 
  scale_x_continuous(breaks = seq(1, ncol(pfm), by = 1)) +
  labs(x = "Position", y = "Frequency", title = "Custom DNA Binding Motif") +
  theme_bw()

pfm <- t(apply(pfm, 1, function(x) x / max(x)))

# Create the sequence logo
seqlogo_obj <- ggseqlogo(pfm, method = "probability")
pfm_normalized <- t(apply(pfm, 1, function(x) x / max(x)))

# Create the sequence logo
seqlogo_obj <- ggseqlogo(pfm_normalized, method = "probability")

# Create a data frame for the variable profile
variable_profile_df <- data.frame(Position = 1:ncol(pfm_normalized), Frequency = apply(pfm_normalized, 2, max))

# Create the variable profile plot
variable_profile_plot <- ggplot(variable_profile_df, aes(x = Position, y = Frequency)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(x = "Position", y = "Normalized Frequency", title = "Custom DNA Binding Motif") +
  theme_void()

# Combine the sequence logo and variable profile using cowplot
library(cowplot)

combined_plot <- plot_grid(seqlogo_obj, variable_profile_plot, nrow = 1, labels = c("A", "B"))

print(combined_plot)





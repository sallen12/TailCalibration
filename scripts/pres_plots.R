################################################################################
##### plot Brehmer and Strokorb mixture distributions

library(texmex)
library(ggplot2)

x <- seq(0, 10, 0.01)
g_x <- dgev(x, mu = 5, sigma = 1, xi = -1/2)
g_x[is.nan(g_x)] <- 0
f_x <- dgev(x, mu = 3, sigma = 1.2, xi = 0)

lambda <- 0.01
f_mix <- lambda * f_x + (1 - lambda) * g_x

df <- data.frame(f = c(g_x, f_mix), x = x, c = rep(c("G", "F"), each = length(x)))
ggplot(df) + geom_ribbon(aes(x = x, ymin = 0, ymax = f, fill = c), alpha = 0.4) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(x[which(g_x == 0)[1]], max(x)),
                     labels = c(expression(x[G]), expression(x[F]))) +
  scale_y_continuous(expand = expansion(c(0, 0.5))) +
  scale_fill_manual(values = c("#F8766d", "#00BFC4"), labels = c(expression(F[lambda]), "G")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ggtitle(expression(paste(lambda, " = 0.01")))
ggsave("plots/score_eg_001.png", width = 3, height = 3.2)


################################################################################
##### plot non-random forecast example

x <- seq(-5, 5, 0.01)
g_x <- dnorm(x)
f_x <- c(dnorm(x[x < 1.3], 0, 2), dnorm(x[x >= 1.3]))

df <- data.frame(f = c(g_x, f_x), x = x, c = rep(c("G", "F"), each = length(x)))
ggplot(df) + geom_ribbon(aes(x = x, ymin = 0, ymax = f, fill = c), alpha = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(c(0, 0.5))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("plots/cal_eg.png", width = 5, height = 3.5)


################################################################################
##### plot unfocused forecast example

x <- seq(-5, 5, 0.01)
g_x <- dnorm(x)
f_xp <- 0.5*dnorm(x, 0, 1) + 0.5*dnorm(x, 1, 1)
f_xm <- 0.5*dnorm(x, 0, 1) + 0.5*dnorm(x, -1, 1)

df <- data.frame(f = c(g_x, f_xp, f_xm), x = x, c = rep(c("G", "F1", "F2"), each = length(x)))
ggplot(df) + geom_ribbon(aes(x = x, ymin = 0, ymax = f, fill = c), alpha = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(c(0, 0.5))) +
  scale_fill_manual(values = c("#F8766d", "#F8766d", "#00BFC4"),
                    labels = c(expression(F["-"]), expression(F["+"]), "G")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("plots/unf_eg.png", width = 5, height = 3.5)


################################################################################
##### plot crps

x <- seq(-5, 5, 0.01)
F_x <- pnorm(x)

y <- 0.5
G_y <- as.numeric(y <= x)

df <- data.frame(x = x, F_x = F_x, G_y = G_y)
ggplot(df) + geom_line(aes(x = x, y = F_x), col = "#00BFC4", linewidth = 1.2) +
  geom_line(aes(x = x, y = G_y), col = "#F8766d", linewidth = 1.2) +
  scale_x_continuous(name = "x", expand = c(0, 0)) +
  scale_y_continuous(name = "F(x)", breaks = c(0, 1),
    limits = c(-0.01, 1.01), expand = c(0, 0)) +
  geom_ribbon(aes(x = x, ymin = pmin(F_x, G_y), ymax = pmax(F_x, G_y)), fill = "#F8766d", alpha = 0.4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("plots/crps_eg.png", width = 5, height = 3.5)


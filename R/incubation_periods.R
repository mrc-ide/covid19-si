# plot incubation period:

x <- seq(0, 30, by = 0.1)

IBM_inc <- dgamma(x = x, 
                  shape = 1.675513,
                  scale = 3.043843)

Lauer_inc <- dgamma(x = x, 
                    shape = 5.807,
                    scale = 0.948)

inc_df <- data.frame(x, IBM_inc, Lauer_inc)

inc_df <- reshape2::melt(inc_df, id = "x")

I <- ggplot(inc_df, aes(x = x + params_offset, y = value, group = variable))+
  geom_line(aes(col = variable), size = 2)+
  scale_color_manual(values = c("red", "blue"), )+
  geom_vline(xintercept = unlist(gamma_shapescale2mucv(shape = 5.807, scale = 0.948)["mu"]) + params_offset,
             col = "blue", size = 2, lty = 2)+
  geom_vline(xintercept = unlist(gamma_shapescale2mucv(shape = 1.675513, scale = 3.0438343)["mu"]) + params_offset,
             col = "red", size = 2, lty = 2)+
  theme_minimal()+
  theme(legend.text=element_text(size=20))

I +
  geom_histogram(
    data = real_data, aes(x = si, y = ..density.., col = NA, group = NA),
    alpha = 0.3,
    binwidth = 1
  )


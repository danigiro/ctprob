library(tsibble)
data <- tsibbledata::pelt
lynxhare <- data |>
  mutate(Year = as.Date(ISOdate(Year, 1,1))) |>
  pivot_longer(-Year) |>
  ggplot(aes(x = Year, y = value/10000, col = name)) + 
  geom_line() +
  labs(y = expression("Number of pelts traded ("*10^4*")")) + 
  theme_bw() + 
  theme(legend.position = "top",
        text = element_text(size = 11),
        legend.margin = margin(b = -10),
        legend.title = element_blank(),
        legend.justification = "right")


ggsave("./analisys/lynxhare.pdf", plot = lynxhare,
       width = 10,
       height = 4)
library(tidyverse)
K <- c(2,1)

extract_cov <- function(comb, prob, meth, fix){
  K <- c(2,1)
  if(missing(comb)){
    name_comb = NA
  }else{
    name_comb = comb
  }
  ctlist <- NULL
  if(missing(fix)){
    files <- sort(list.files(file.path("./ProbReco", prob, meth), full.names = TRUE))
  }else{
    files <- sort(list.files(file.path("./ProbReco", prob, meth), full.names = TRUE))[fix]
  }
  for(i in 1:length(files)){
    load(files[i])
    if(missing(comb)){
      res <- t(do.call(rbind, listBase[paste0("k", sort(K, decreasing = TRUE))]))
    }else{
      res <- t(do.call(rbind, listFree[[comb]][paste0("k", sort(K, decreasing = TRUE))]))
    }
    N <- NCOL(res) / sum(K)
    DN <- FoReco:::Dmat(h = N, m = K, n = NROW(res))
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
    ctlist[[i]] <- E
  }
  cov <- cov(Reduce(rbind, ctlist))
  rbind(as_tibble(cov, rownames = "row") |>
    pivot_longer(-row, names_to = "col") |>
    mutate_if(is.character, function(x) readr::parse_number(x)) |>
    add_column(comb = name_comb,
               meth = meth,
               prob = prob,
               type = "cov"),
    as_tibble(cov2cor(cov), rownames = "row") |>
      pivot_longer(-row, names_to = "col") |>
      mutate_if(is.character, function(x) readr::parse_number(x)) |>
      add_column(comb = name_comb,
                 meth = meth,
                 prob = prob,
                 type = "cor"))
}
ctjb_base <- extract_cov(prob = "ctjb", meth = "base", fix = 1)
ctsam_base <- extract_cov(prob = "ctsam", meth = "base", fix = 1)
ctsamh_base <- extract_cov(prob = "ctsamh", meth = "base", fix = 1)
hbsam_base <- extract_cov(prob = "hbsam", meth = "base", fix = 1)
hbsamh_base <- extract_cov(prob = "hbsamh", meth = "base", fix = 1)

ctjb_oct_sam <- extract_cov(comb = "sam", prob = "ctjb", meth = "oct", fix = 1)
ctsam_oct_sam <- extract_cov(comb = "sam", prob = "ctsam", meth = "oct", fix = 1)
ctsamh_oct_sam <- extract_cov(comb = "sam", prob = "ctsamh", meth = "oct", fix = 1)
ctjb_octh_sam <- extract_cov(comb = "sam", prob = "ctjb", meth = "octh", fix = 1)
ctsam_octh_sam <- extract_cov(comb = "sam", prob = "ctsam", meth = "octh", fix = 1)
ctsamh_octh_sam <- extract_cov(comb = "sam", prob = "ctsamh", meth = "octh", fix = 1)
cov_sam_figure <- rbind(ctjb_oct_sam, ctsam_oct_sam, hbsam_base, hbsamh_base,
      ctsamh_oct_sam, ctjb_octh_sam, ctsam_octh_sam, ctsamh_octh_sam,
      ctsam_base, ctjb_base, ctsamh_base) |>
  filter(type =="cor") |>
  ggplot(aes(x = row, y = col, fill = value)) + 
  geom_tile(col = "black") +
  scale_fill_gradient2() +
  scale_y_reverse() + 
  geom_text(aes(label = round(value, 2)))+
  theme_void() +
  facet_wrap(vars(prob, meth, comb), nrow = 3, dir = "v")+
  coord_fixed()

save(ctjb_base, ctsamh_base, ctsam_base, file = "analisys/base_cov_plot.RData")

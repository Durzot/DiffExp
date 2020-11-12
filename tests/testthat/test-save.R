test_that("update table works", {
  tab_1 <- data.frame(meta="meta", b=c(1,2,3))
  tab_2 <- data.frame(meta="meta", b=c(4,5,6))

  tab <- update_table(tab_1, tab_2, tags=list(meta="meta"))
  expect_true(all.equal(tab, tab_2))

  tab_1 <- data.frame(meta="meta_1", b=c(1,2,3))
  tab_2 <- data.frame(meta="meta_2", b=c(4,5,6))

  tab <- update_table(tab_1, tab_2, tags=list(meta="meta"))
  expect_true(all.equal(tab, rbind(tab_1,tab_2)))
})

test_that("save update table works", {
  tab_1 <- data.frame(meta="meta", b=c(1,2,3))
  tab_2 <- data.frame(meta="meta", b=c(4,5,6))
  file <- "./out/test_save_update.txt"

  save_update_table(file, tab_1, tags=list(meta="meta"), verbose=F)
  save_update_table(file, tab_2, tags=list(meta="meta"), verbose=F)
  tab <- read.table(file, sep="\t", header=T)
  expect_true(all.equal(tab, tab_2))

  save_update_table(file, tab_2, tags=list(meta="meta"), verbose=F)
  tab <- read.table(file, sep="\t", header=T)
  expect_true(all.equal(tab, tab_2))
})

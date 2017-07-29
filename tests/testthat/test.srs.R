context('successor representations')

D <- simulateFixed(20, 10, 2)

expect_identical(comp.SRSA(D$D, 6, 0.2, 0.2),
                 matrix(srs(D$D, 6, 0.2, 0.2, normalize=TRUE),
                        ncol=36)
                 )



expect_equal(CostSRSA(c(0.2, 0.5), D$D, 6, D$ppt$score, 2, TRUE),
             costPCA(c(0.2, 0.5), D$D, 6, D$ppt$score, 2, NULL, TRUE, FALSE)
                 )

system.time(CostSRSA(c(0.2, 0.5), D$D, 6, D$ppt$score, 2, TRUE))
system.time(costPCA(c(0.2, 0.5), D$D, 6, D$ppt$score, 2, NULL, TRUE, FALSE))



microbenchmark(
  srsaPCA(D$D, 6, D$ppt$score, runif(1), runif(1), 2, NULL),
  DoSRSA(D$D, 6, D$ppt$score, runif(1), runif(1), 2)
)
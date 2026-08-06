source("r0numerical2specIntro.r", chdir=FALSE)
source("r0numerical2specIntroOldDenom.r", chdir=FALSE)

BB <- 0.15
cross <- 0.0005

args <- list(
  bm1h=BB, bhm1=BB, bm1p1=BB, bm2h=BB, bhm2=BB, bm2p1=BB,
  bp1m2=BB, bp1m1=BB,
  bm2p2=BB, bp2m1=BB, bm1p2=BB, bp2m2=BB,
  rm1h=0.5, rm1p1=cross, rm2h=cross, rm2p1=0.5,
  rm2p2=0.5, rm1p2=0.5,
  Vm=1/7, Vm2=1/7,
  gh=1/4, mh=1/(60*365), gp1=1/4, mp1=1/(15*365),
  mp2=1/(30*365), gp2=1/5,
  Nh=1000, Np1=1000, Np2=0, Nm1=25000, Nm2=25000,
  NN=2000, intro=0
)

R0_new <- do.call(r0new, args)
R0_old <- do.call(r0old, args)

Nh_vec <- c(args$Nh, args$Np1)
r_v1 <- c(args$rm1h, args$rm1p1)
r_v2 <- c(args$rm2h, args$rm2p1)

rescale <- function(r_vec, Nh_vec) {
  c_j <- sum(Nh_vec) / sum(r_vec * Nh_vec)
  r_vec * c_j
}

tilde_v1 <- rescale(r_v1, Nh_vec)
tilde_v2 <- rescale(r_v2, Nh_vec)

args_eq <- args
args_eq$rm1h <- tilde_v1[1]
args_eq$rm1p1 <- tilde_v1[2]
args_eq$rm2h <- tilde_v2[1]
args_eq$rm2p1 <- tilde_v2[2]

R0_old_eq <- do.call(r0old, args_eq)

message("=== Equal host populations (Nh=Np1=1000) ===")
message("R0_new (weighted, original r): ", R0_new)
message("R0_old (unweighted, original r): ", R0_old)
message("R0_old_eq (unweighted, rescaled r): ", R0_old_eq)
message("Diff R0_new vs R0_old_eq: ", abs(R0_new - R0_old_eq))
message("c_1 = ", sum(Nh_vec)/sum(r_v1*Nh_vec))
message("c_2 = ", sum(Nh_vec)/sum(r_v2*Nh_vec))

# Test 2: unequal hosts
args2 <- args
args2$Nh <- 2500; args2$Np1 <- 500; args2$NN <- 3000

R0_new2 <- do.call(r0new, args2)
R0_old2 <- do.call(r0old, args2)

Nh_vec2 <- c(2500, 500)
tilde2_v1 <- rescale(r_v1, Nh_vec2)
tilde2_v2 <- rescale(r_v2, Nh_vec2)

args2_eq <- args2
args2_eq$rm1h <- tilde2_v1[1]; args2_eq$rm1p1 <- tilde2_v1[2]
args2_eq$rm2h <- tilde2_v2[1]; args2_eq$rm2p1 <- tilde2_v2[2]

R0_old2_eq <- do.call(r0old, args2_eq)

message("")
message("=== Unequal host populations (Nh=2500, Np1=500) ===")
message("R0_new: ", R0_new2)
message("R0_old: ", R0_old2)
message("R0_old_eq: ", R0_old2_eq)
message("Diff R0_new vs R0_old_eq: ", abs(R0_new2 - R0_old2_eq))

# Test 3: very asymmetric 
args3 <- args
args3$Nh <- 4000; args3$Np1 <- 200; args3$NN <- 4200

R0_new3 <- do.call(r0new, args3)

Nh_vec3 <- c(4000, 200)
tilde3_v1 <- rescale(r_v1, Nh_vec3)
tilde3_v2 <- rescale(r_v2, Nh_vec3)

args3_eq <- args3
args3_eq$rm1h <- tilde3_v1[1]; args3_eq$rm1p1 <- tilde3_v1[2]
args3_eq$rm2h <- tilde3_v2[1]; args3_eq$rm2p1 <- tilde3_v2[2]

R0_old3_eq <- do.call(r0old, args3_eq)

message("")
message("=== Very asymmetric (Nh=4000, Np1=200) ===")
message("R0_new: ", R0_new3)
message("R0_old_eq: ", R0_old3_eq)
message("Diff: ", abs(R0_new3 - R0_old3_eq))

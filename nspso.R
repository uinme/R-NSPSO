currentDirectory <- (function() {
   path <- (function() attr(body(sys.function(-1)), "srcfile"))()$filename
   dirname(path)
})()
setwd(currentDirectory)
source("lib/non_dominated.R")
source("lib/crowding_distance_sort.R")

if (!file.exists("image")) {
	dir.create("image")
}

AN <- 100       # Number of particle
FN <- 2         # Number of objective
DVN <- 3        # Number of design variable
DIM <- DVN # DVM * EVN
C1 <- 2.05
C2 <- 2.05
LW <- 0.1
UW <- 0.9
DVLT <- rep(0 ,DIM)
DVUT <- rep(1, DIM)
Itmax <- 100

vMax <- (DVUT - DVLT) / 2
vMin <- -vMax

# Objective function
# ZDT1 function
fn <- function(x) {
	n <- length(x)
	g <- 1.0 + (9 / (n -1)) * sum(x[2:n])

	f1 <- x[1]
	f2 <- g * (1.0 - sqrt(f1 / g))

	return(c(f1, f2))
}

# Definition of variables
# Population
popFit   <- matrix(0, nrow=AN, ncol=FN)
popVar   <- matrix(0, nrow=AN, ncol=DIM)
pbestFit <- matrix(0, nrow=AN, ncol=FN)
pbestVar <- matrix(0, nrow=AN, ncol=DIM)
gbestFit <- matrix(0, nrow=AN, ncol=FN)
gbestVar <- matrix(0, nrow=AN, ncol=DIM)
velocity <- matrix(0, nrow=AN, ncol=DIM)
# Temporary non-dominated
ndomFit <- matrix(0, nrow=2*AN, ncol=FN)
ndomVar <- matrix(0, nrow=2*AN, ncol=DIM)
ndomCnt <- 0
# Initialize progress bar
pb <- txtProgressBar(min=1, max=Itmax, style=3)

# Initialize the population
#   Randomize the variable of each particles.
for (i in 1:DIM) {
	popVar[,i] <- runif(AN, DVLT[i], DVUT[i])
}
popFit[1:AN, ] <- t(apply(popVar, 1, fn))
pbestFit <- popFit
pbestVar <- popVar

Itn <- 0
while (Itn < Itmax) {
	Itn <- Itn + 1

	# Evaluation
	popFit[1:AN, ] <- t(apply(popVar, 1, fn))

	# Non-dominate
	unionFit <- rbind(popFit, pbestFit)
	unionVar <- rbind(popVar, pbestVar)
	ixndom   <- non_dominated(unionFit, 0)
	ndomCnt  <- length(ixndom)

	# Sort the non-dominated particles on basis of crowding distance
	ndomFit[1:ndomCnt, ] <- unionFit[ixndom, ]
	ndomVar[1:ndomCnt, ] <- unionVar[ixndom, ]
	ixcrow  <- crowding_distance_sort(ndomFit[1:ndomCnt, ], "min")
	ndomFit[1:ndomCnt, ] <- ndomFit[ixcrow, ]
	ndomVar[1:ndomCnt, ] <- ndomVar[ixcrow, ]

	# Render the objective space and the variable space
	fname <- paste("image/image", formatC(Itn, width=3, flag="0"), ".png", sep="")
	png(filename=fname, width=800, height=400, pointsize=12, bg="white")
	par(mfrow=c(1, 2))
	# The objective space
	plot(popFit[ ,1], popFit[ ,2], xlim=c(0, 2), ylim=c(0, 2), col="gray", pch=1,
			 main="Objective Space")
	points(ndomFit[1:ndomCnt ,1], ndomFit[1:ndomCnt ,2], col="cyan4", pch=16)
	points(gbestFit[ ,1], gbestFit[ ,2], col="red", pch=1)
	points(pbestFit[ ,1], pbestFit[ ,2], col="blue", pch=1)
	legend("topright", legend=c("Particle", "non-dominated", "gbest", "pbest"),
				 col=c("gray", "cyan4", "red", "blue"), pch=c(1, 16, 1, 1))

	# The variable space
	plot(popVar[ ,1], popVar[ ,2], xlim=c(0, 2), ylim=c(0, 2), col="gray", pch=1,
			 main="Variable Space")
	points(ndomVar[1:ndomCnt ,1], ndomVar[1:ndomCnt ,2], col="cyan4", pch=16)
	points(gbestVar[ ,1], gbestVar[ ,2], col="red", pch=1)
	points(pbestVar[ ,1], pbestVar[ ,2], col="blue", pch=1)
	legend("topright", legend=c("Particle", "non-dominated", "gbest", "pbest"),
	       col=c("gray", "cyan4", "red", "blue"), pch=c(1, 16, 1, 1))
	dev.off()
	
	# Update the gbest
	top <- ceiling(ndomCnt * 0.05)
	ixgbest <- sample(1:top, AN, replace=TRUE)
	gbestFit[1:AN, ] <- ndomFit[ixgbest, ]
	gbestVar[1:AN, ] <- ndomVar[ixgbest, ]

	# Update the pbest
	if (ndomCnt > AN) {
		pbestCnt <- AN
		pbestFit[1:AN, ] <- unionFit[ixndom[1:AN], ]
		pbestVar[1:AN, ] <- unionVar[ixndom[1:AN], ]
	}
	else {
		pbestCnt <- ndomCnt
		pbestFit[1:ndomCnt, ] <- unionFit[ixndom, ]
		pbestVar[1:ndomCnt, ] <- unionVar[ixndom, ]
		restFit <- unionFit[-ixndom, ]
		restVar <- unionVar[-ixndom, ]
		while (pbestCnt < AN) {
			ixndomr  <- non_dominated(restFit, 0)
			ndomrCnt <- length(ixndomr)
			needCnt  <- AN - pbestCnt
			if (needCnt < ndomrCnt) {
				pbestFit[(pbestCnt+1):AN, ] <- restFit[ixndomr[1:needCnt], ]
				pbestVar[(pbestCnt+1):AN, ] <- restVar[ixndomr[1:needCnt], ]
				pbestCnt <- AN
			}
			else {
				pbestFit[(pbestCnt+1):(pbestCnt+ndomrCnt), ] <- restFit[ixndomr, ]
				pbestVar[(pbestCnt+1):(pbestCnt+ndomrCnt), ] <- restVar[ixndomr, ]
				pbestCnt <- pbestCnt + ndomrCnt
			}
			restFit <- restFit[-ixndomr, ]
			restVar <- restVar[-ixndomr, ]
		}
	}

	# Update the velocity and position of each particle
	W  <- UW - ((UW - LW) / Itmax) * Itn
	R1 <- matrix(runif(rep(AN, DIM), 0, 1), nrow=AN, ncol=DIM)
	R2 <- matrix(runif(rep(AN, DIM), 0, 1), nrow=AN, ncol=DIM)
	velocity[1:AN, ] <- 0.7298438 * (W * velocity[1:AN, ] +
																	 C1 * R1 * (pbestVar[1:AN, ] - popVar[1:AN, ]) +
																	 C2 * R2 * (gbestVar[1:AN, ] - popVar[1:AN, ]))

	# Limit the velocity
	for (i in 1:DIM) {
		ixmax <- vMax[i] < velocity[1:AN, i]
		ixmin <- vMin[i] > velocity[1:AN, i]
		velocity[ixmax, i] <- vMax[i]
		velocity[ixmin, i] <- vMin[i]
	}

	# Update the position of each particle
	popVar[1:AN, ] <- popVar[1:AN, ] + velocity[1:AN, ]

	# Bound (Reflection)
	for (i in 1:DIM) {
		ixmax <- DVUT[i] < popVar[1:AN, i]
		ixmin <- DVLT[i] > popVar[1:AN, i]
		popVar[ixmax, i] <- popVar[ixmax, i] - velocity[ixmax, i]
		popVar[ixmin, i] <- popVar[ixmin, i] - velocity[ixmin, i]
	}


	setTxtProgressBar(pb, Itn)
}


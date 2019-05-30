LOG <- function (par, unif, structure) 
{
  tt <- NULL
  if (length(unique(unif)) != length(unif)) 
    stop("The 'unif' argument must be composed of different values")
  if (par > 1 || par < 0) 
    stop("Wrong 'param' input")
  if (is.null(structure)) {
    t <- new("Log_Child", parameter = as.character(par), 
             arg = unif, dimension = length(unif), name = "Logarithmic distribution", 
             type = "Child", obj = "Log")
  }
  else {
    if (class(structure) != "list") 
      stop("The argument 'structure' must be a list")
    if (is.null(unif)) 
      t <- new("Log_Mother", parameter = as.character(par), 
               dimension = length(structure), structure = structure, 
               arg = 0, name = "Logarithmic distribution", type = "Mother", 
               obj = "Log")
    else t <- new("Log_Mother", parameter = as.character(par), 
                  dimension = length(structure) + length(unif), structure = structure, 
                  arg = unif, name = "Logarithmic distribution", type = "Mother", 
                  obj = "Log")
  }
  if (t@type == "Mother") {
    t@Param <- "gamma"
    t@Laplace <- "log(1 - (gamma) * exp(-(z))) / log(1 - (gamma))"
    t@LaplaceInv <- "-log((1 - (1 - (gamma))^(z)) / (gamma))"
    t@PGF <- "log(1 - (gamma)*(z)) / log(1 - (gamma))"
    t@PGFInv <- "((1 - (1 - (gamma))^(z))/(gamma))"
    t@Der <- function(tt, k, type) {
      if (type == "PGF") {
        if (k >= 1) {
          ini <- paste("(", -factorial(k - 1), ") / log(1 - gamma) * (gamma / (1 - gamma * (z)))^(", 
                       k, ")", sep = "")
          stringr::str_replace_all(ini, "z", tt)
        }
        else t@PGF
      }
      else if (type == "PGFInv") {
        if (k == 1) {
          ini <- "-log(1 - gamma) / gamma * (1 - gamma)^(z)"
          stringr::str_replace_all(ini, "z", tt)
        }
        else if (k == 0) 
          t@PGFInv
      }
    }
    t@FUN <- function(type) {
      if (type == "PGF") 
        function(tt, gamma) log(1 - (gamma) * (tt))/log(1 - 
                                                          (gamma))
      else if (type == "PGFInv") 
        function(tt, gamma) ((1 - (1 - (gamma))^(tt))/(gamma))
      else if (type == "PGF.Der") {
        function(tt, gamma, k) -factorial(k - 1)/log(1 - 
                                                       gamma) * (gamma/(1 - gamma * (tt)))^(k)
      }
      else if (type == "PGFInv.Der") {
        function(tt, gamma) -log(1 - gamma)/gamma * (1 - 
                                                       gamma)^(tt)
      }
    }
  }
  else {
    t@Param <- "alpha"
    t@Laplace <- "log(1 - (alpha) * exp(-(z))) / log(1 - (alpha))"
    t@LaplaceInv <- "-log((1 - (1 - (alpha))^(z)) / (alpha))"
    t@PGF <- "log(1 - (alpha)*(z)) / log(1 - (alpha))"
    t@PGFInv <- "((1 - (1 - (alpha))^(z))/(alpha))"
    t@Der <- function(tt, k, type) {
      if (type == "PGF") {
        if (k >= 1) {
          ini <- paste("(", -factorial(k - 1), ") / log(1 - alpha) * (alpha / (1 - alpha * (z)))^(", 
                       k, ")", sep = "")
          stringr::str_replace_all(ini, "z", tt)
        }
        else stringr::str_replace_all(t@PGF, "z", tt)
      }
      else if (type == "PGFInv") {
        if (k == 1) {
          ini <- "-log(1 - alpha) / alpha * (1 - alpha)^(z)"
          stringr::str_replace_all(ini, "z", tt)
        }
        else if (k == 0) 
          stringr::str_replace_all(t@PGFInv, "z", tt)
      }
      else if (type == "Laplace") {
        if (k >= 1) {
          res <- numeric(k)
          for (r in 1:k) {
            ini <- t@Der("exp(-(z))", r, "PGF")
            ini <- paste("(", ini, ") * (exp(-", r, " * (z)) * (-1)^(", 
                         k, "))", sep = "")
            input <- (-1)^(r - 1)/factorial(1)/factorial(r - 
                                                           1) * 1^k
            if (r > 1) {
              for (s in 2:r) {
                input <- input + ((-1)^(r - s)/factorial(s)/factorial(r - 
                                                                        s) * s^k)
              }
            }
            res[r] <- paste("(", ini, ") * (", input, 
                            ")", sep = "")
          }
          res <- paste("(", res, ")", sep = "", collapse = " + ")
          stringr::str_replace_all(res, "z", tt)
        }
        else {
          stringr::str_replace_all(t@Laplace, "z", tt)
        }
      }
      else if (type == "LaplaceInv") {
        ini <- paste("-(", t@Der(tt, 1, "PGFInv"), ") / (", 
                     stringr::str_replace_all(t@PGFInv, "z", tt), 
                     ")", sep = "")
        ini
      }
    }
    t@FUN <- function(type) {
      if (type == "PGF") 
        function(tt, gamma) log(1 - (gamma) * (tt))/log(1 - 
                                                          (gamma))
      else if (type == "PGFInv") 
        function(tt, gamma) ((1 - (1 - (gamma))^(tt))/(gamma))
      else if (type == "Laplace") 
        log(1 - (gamma) * exp(-(tt)))/log(1 - (gamma))
      else if (type == "LaplaceInv") 
        -log((1 - (1 - (gamma))^(tt))/(gamma))
      else if (type == "PGF.Der") {
        function(tt, gamma, k) -factorial(k - 1)/log(1 - 
                                                       gamma) * (gamma/(1 - gamma * (tt)))^(k)
      }
      else if (type == "PGFInv.Der") {
        function(tt, gamma) -log(1 - gamma)/gamma * (1 - 
                                                       gamma)^(tt)
      }
    }
  }
  t@simul <- function(z, gamma) copula::rlog(z, gamma)
  t@theta <- vector("numeric")
  t@cop <- function(gamma, dim) Frank(gamma, dim)
  t
}


LOG(0.5,1,NULL)

rComSim <- function (n, structure) 
{
  e1 <- new.env()
  e1$res <- list()
  gen <- GeneticCodes(structure)
  e1$M0 <- structure@simul(n, as.numeric(structure@parameter))
  for (i in 1:length(gen)) {
    if (length(gen[[i]]) == 2) {
      str2 <- Node(gen[[i]], structure)
      M.prec <- eval(parse(text = paste("e1$M", paste(gen[[i]][1], 
                                                      collapse = ""), sep = "")))
      R <- matrix(rexp(length(str2@arg) * n, 1), ncol = length(str2@arg), 
                  nrow = n)
      if (gen[[i]][length(gen[[i]])] == 0) {
        Theta <- M.prec
        ini <- stringr::str_replace_all(structure@Laplace, 
                                        structure@Param, structure@parameter)
        ff <- function(z) eval(parse(text = ini))
        e1$res[[i]] <- ff(R/Theta)
      }
      else {
        Theta <- matrix(rep(vapply(1:length(M.prec), 
                                   function(t) sum(str2@simul(M.prec[t], as.numeric(str2@parameter))), 
                                   0), length(str2@arg)), ncol = length(str2@arg), 
                        nrow = n)
        ini <- stringr::str_replace_all(structure@PGF, 
                                        structure@Param, structure@parameter)
        ini <- stringr::str_replace_all(ini, "z", stringr::str_replace_all(str2@Laplace, 
                                                                           str2@Param, str2@parameter))
        ff <- function(z) eval(parse(text = ini))
        e1$res[[i]] <- ff(R/Theta)
      }
    }
    else if (length(gen[[i]]) > 2) {
      Lap <- stringr::str_replace_all(structure@PGF, structure@Param, 
                                      structure@parameter)
      for (j in 2:(length(gen[[i]]) - 1)) {
        str2 <- Node(gen[[i]][1:j], structure)
        if (gen[[i]][length(gen[[i]])] != 0) {
          ini <- stringr::str_replace_all(str2@PGF, str2@Param, 
                                          str2@parameter)
          Lap <- stringr::str_replace_all(Lap, "z", ini)
        }
        else {
          if (j == length(gen[[i]]) - 1) {
            ini <- stringr::str_replace_all(str2@Laplace, 
                                            str2@Param, str2@parameter)
            Lap <- stringr::str_replace_all(Lap, "z", 
                                            ini)
          }
          else {
            ini <- stringr::str_replace_all(str2@PGF, 
                                            str2@Param, str2@parameter)
            Lap <- stringr::str_replace_all(Lap, "z", 
                                            ini)
          }
        }
        variable0 <- paste("M", paste(gen[[i]][1:(j - 
                                                    1)], collapse = ""), sep = "")
        variable1 <- paste("M", paste(gen[[i]][1:j], 
                                      collapse = ""), sep = "")
        if (!exists(variable1, envir = e1)) {
          eval(parse(text = paste("e1$", variable1, " <- vapply(1:length(", 
                                  paste("e1$", variable0, sep = ""), "),\n                                  function(t) sum(str2@simul(", 
                                  paste("e1$", variable0, "[t],", sep = ""), 
                                  "as.numeric(str2@parameter))), 0)", sep = "")))
        }
      }
      str2 <- Node(gen[[i]], structure)
      M.prec <- eval(parse(text = paste("e1$M", paste(gen[[i]][1:(length(gen[[i]]) - 
                                                                    1)], collapse = ""), sep = "")))
      if (gen[[i]][length(gen[[i]])] != 0) {
        ini <- stringr::str_replace_all(str2@Laplace, 
                                        str2@Param, str2@parameter)
        Lap <- stringr::str_replace_all(Lap, "z", ini)
        Theta <- matrix(rep(vapply(1:length(M.prec), 
                                   function(t) sum(str2@simul(M.prec[t], as.numeric(str2@parameter))), 
                                   0), length(str2@arg)), ncol = length(str2@arg), 
                        nrow = n)
      }
      else Theta <- M.prec
      R <- matrix(rexp(length(str2@arg) * n, 1), ncol = length(str2@arg), 
                  nrow = n)
      ff <- function(z) eval(parse(text = Lap))
      e1$res[[i]] <- ff(R/Theta)
    }
  }
  do.call(cbind, e1$res)
}

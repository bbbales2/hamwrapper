#' Create a hamiltonian given a Stan model and data for that model. Can optionally specify a metric
#' 
#' The hamiltonian is just a list with a bunch of named elements. These elements are:
#' 
#' \describe{
#'   \item{M}{Euclidean metric} 
#'   \item{data}{A copy of the data passed to the model}
#'   \item{env}{The R environment where the compiled code is stored}
#'   \item{U}{Function that takes one argument, the position, and returns the potential energy}
#'   \item{H}{Function that takes two arguments, the position and momentum, and returns the value of the Hamiltonian}
#'   \item{gradU}{Function that takes one argument, the position, and returns the gradient of the potential energy}
#'   \item{hessU}{Function that takes one argument, the position, and returns the hessian of the potential energy}
#'   \item{hessUVecProd}{Function that takes two arguments, the position and a vector, and returns the product of the hessian
#'      of the potential energy and the vector}
#'   \item{sampleMomentum}{Function that takes no arguments, but returns a random momentum sample suitable for HMC}
#' }
#'
#' @param modelFile Stan model
#' @param data Data for model
#' @param M Euclidean metric for hamiltonian
#'
#' @return Hamiltonian list
#' @export
#'
createHamiltonianSystem <- function(modelFile, data, M = NULL) {
  env = new.env()
  
  data_file = tempfile()
  if(length(data) == 0) {
    file.create(data_file)
  } else {
    rstan::stan_rdump(names(data), file = data_file, envir = list2env(data))
  }
  
  # create linear_regression_model.cpp from Stan file
  if(!file.exists(modelFile)) {
    stop("model '", modelFile,"' does not exist")
  }
  code = rstan::stanc(modelFile, model_name = 'log_density', obfuscate_model_name = FALSE)
  
  # compile with RCPP
  helper_file = system.file("extdata", "gradient_helper.cpp", package = "hamwrapper")
  helper_code = readChar(helper_file, file.info(helper_file)$size)
  model_code_file = tempfile(fileext = ".cpp")
  write(paste(code$cppcode, helper_code, sep = "\n"), model_code_file)
  Rcpp::sourceCpp(model_code_file, env = env)
  env$set_data(data_file)
  file.remove(data_file)
  file.remove(model_code_file)
  
  if(is.null(M)) {
    M = diag(env$get_num_params())
  }
  
  L = chol(M)
  
  return(list(M = M,
             data = data,
             env = env,
             U = function(q) { -env$jacobian(q)$u },
             H = function(q, p) { 0.5 * sum(p * solve(M, p)) - env$jacobian(q)$u },
             gradU = function(q) { -env$jacobian(q)$jac },
             hessU = function(q) { -env$hessian(q)$hess },
             hessUVecProd = function(q, vec) { -env$hessian_vector(q, vec)$hessv },
             sampleMomentum = function() { L %*% rnorm(nrow(L)) }))
}

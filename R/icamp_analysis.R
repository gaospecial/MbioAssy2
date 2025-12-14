#' Run iCAMP Analysis
#'
#' @param comm Community abundance matrix (rows are samples, columns are OTUs).
#' @param tree Phylogenetic tree object (phylo).
#' @param treat Group information (data frame).
#' @param env Environmental variables (data frame, optional).
#' @param clas Taxonomy information (data frame).
#' @param prefix Prefix for output files.
#' @param save_dir Directory to save outputs and intermediate files.
#' @param rand Randomization times, default 1000.
#' @param nworker Number of threads, default 4.
#' @param memory_G Memory limit in GB, default 50.
#' @param bin_size_limit Minimum bin size, default 24.
#' @return List containing iCAMP results (icres and icbin).
#' @importFrom iCAMP match.name pdist.big dniche taxa.binphy.big ps.bin icamp.big icamp.bins
#' @importFrom utils write.csv
#' @importFrom ape read.tree
#' @export
run_icamp <- function(comm, tree, treat, env = NULL, clas, prefix = "iCAMP_Test", save_dir = tempdir(), 
                      rand = 1000, nworker = 4, memory_G = 50, bin_size_limit = 24) {
  
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # Preprocessing: Match IDs
  # If env is NULL, match without it
  if (!is.null(env)) {
    sampid.check <- iCAMP::match.name(rn.list = list(comm = comm, treat = treat, env = env))
    env <- sampid.check$env
  } else {
    sampid.check <- iCAMP::match.name(rn.list = list(comm = comm, treat = treat))
  }
  treat <- sampid.check$treat
  comm <- sampid.check$comm
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  
  spid.check <- iCAMP::match.name(
    cn.list = list(comm = comm),
    rn.list = list(clas = clas),
    tree.list = list(tree = tree)
  )
  comm <- spid.check$comm
  clas <- spid.check$clas
  tree <- spid.check$tree
  
  # 1. Phylogenetic Distance Matrix
  pd_wd <- file.path(save_dir, "pd")
  if (!dir.exists(pd_wd)) dir.create(pd_wd)
  
  pd.big <- iCAMP::pdist.big(
    tree = tree,
    wd = pd_wd,
    nworker = nworker,
    memory.G = memory_G
  )
  
  # 2. iCAMP Analysis
  # Note: The original script does niche diff and phylo signal before this.
  # For the wrapper, we go straight to iCAMP main analysis to keep it simple, 
  # assuming the user might check binning separately if needed.
  # Or we could include it if env is provided, but icamp.big doesn't strictly require it 
  # unless we are optimizing binning. Here we use fixed bin.size.limit.
  
  icres <- iCAMP::icamp.big(
    comm = comm,
    pd.desc = pd.big$pd.file,
    pd.spname = pd.big$tip.label,
    pd.wd = pd.big$pd.wd,
    rand = rand,
    tree = tree,
    prefix = prefix,
    ds = 0.2, # Default from script
    pd.cut = NA,
    sp.check = TRUE,
    phylo.rand.scale = "within.bin",
    taxa.rand.scale = "across.all",
    phylo.metric = "bMPD",
    sig.index = "Confidence",
    bin.size.limit = bin_size_limit,
    nworker = nworker,
    memory.G = memory_G,
    rtree.save = FALSE,
    detail.save = TRUE,
    qp.save = FALSE,
    detail.null = FALSE,
    ignore.zero = TRUE,
    output.wd = save_dir,
    correct.special = TRUE,
    unit.sum = rowSums(comm),
    special.method = "depend",
    ses.cut = 1.96,
    rc.cut = 0.95,
    conf.cut = 0.975,
    omit.option = "no",
    meta.ab = NULL
  )
  
  # 3. Bin Level Statistics
  icbin <- iCAMP::icamp.bins(
    icamp.detail = icres$detail,
    treat = treat,
    clas = clas,
    silent = FALSE,
    boot = TRUE,
    rand.time = rand,
    between.group = TRUE
  )
  
  # Save main summaries to files in save_dir as in original script
  utils::write.csv(icbin$Pt, file = file.path(save_dir, paste0(prefix, ".ProcessImportance_EachGroup.csv")), row.names = FALSE)
  utils::write.csv(icbin$Ptk, file = file.path(save_dir, paste0(prefix, ".ProcessImportance_EachBin_EachGroup.csv")), row.names = FALSE)
  
  return(list(icres = icres, icbin = icbin))
}

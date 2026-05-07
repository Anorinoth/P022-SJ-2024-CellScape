#' Clean Up Orphaned Parallel Workers
#' Author: Dr. Michael Zaiken
#'
#' When R scripts crash, BiocParallel SnowParam workers can persist as zombie processes,
#' consuming RAM. This function finds and terminates orphaned worker processes.
#'
#' @param kill_all Logical, if TRUE kills ALL R worker processes (use with caution).
#'   If FALSE, only attempts to stop currently registered BiocParallel clusters.
#'   Default: FALSE
#' @param verbose Logical, whether to print detailed information. Default: TRUE
#'
#' @return List with:
#'   \item{stopped_clusters}{Number of BiocParallel clusters stopped}
#'   \item{killed_processes}{Number of system processes killed (if kill_all=TRUE)}
#'   \item{ram_freed}{Estimated RAM freed in GB}
#'
#' @examples
#' # Safe cleanup (only stops registered clusters)
#' cleanup_parallel_workers()
#'
#' # Aggressive cleanup (kills all R worker processes)
#' cleanup_parallel_workers(kill_all = TRUE)
#'
#' @export
cleanup_parallel_workers <- function(kill_all = FALSE, verbose = TRUE) {
  
  if (verbose) {
    cat("\n=== Cleaning Up Parallel Workers ===\n")
  }
  
  stopped_clusters <- 0
  killed_processes <- 0
  ram_before <- sum(gc()[,2]) * 1048576
  
  # ============================================================================
  # Step 1: Stop currently registered BiocParallel clusters
  # ============================================================================
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    if (verbose) cat("\n1. Checking registered BiocParallel backends...\n")
    
    tryCatch({
      bpp <- BiocParallel::bpparam()
      bpp_class <- class(bpp)[1]
      
      if (verbose) cat("   Current backend:", bpp_class, "\n")
      
      if (inherits(bpp, "SnowParam")) {
        n_workers <- BiocParallel::bpnworkers(bpp)
        if (verbose) cat("   Found SnowParam with", n_workers, "workers\n")
        
        tryCatch({
          parallel::stopCluster(BiocParallel::bpbackend(bpp))
          stopped_clusters <- stopped_clusters + 1
          if (verbose) cat("   ✓ Stopped SnowParam cluster\n")
        }, error = function(e) {
          if (verbose) cat("   ⚠️  Could not stop cluster:", e$message, "\n")
        })
        
        # Reset to SerialParam to prevent reuse of dead cluster
        BiocParallel::register(BiocParallel::SerialParam())
        if (verbose) cat("   ✓ Reset to SerialParam\n")
      } else {
        if (verbose) cat("   No SnowParam cluster found (backend is", bpp_class, ")\n")
      }
    }, error = function(e) {
      if (verbose) cat("   Error checking BiocParallel:", e$message, "\n")
    })
  }
  
  # ============================================================================
  # Step 2: Kill orphaned R worker processes (optional, aggressive)
  # ============================================================================
  if (kill_all) {
    if (verbose) cat("\n2. Killing orphaned R worker processes (kill_all=TRUE)...\n")
    
    # Get current process ID (to avoid killing ourselves)
    current_pid <- Sys.getpid()
    
    # Find all R processes
    ps_cmd <- paste0("ps aux | grep -E 'R|Rscript' | grep -v grep | grep -v '", current_pid, "'")
    ps_output <- system(ps_cmd, intern = TRUE)
    
    if (length(ps_output) > 0) {
      if (verbose) {
        cat("   Found", length(ps_output), "R processes:\n")
        for (line in ps_output) {
          cat("   ", line, "\n")
        }
      }
      
      # Extract PIDs and check for worker processes
      for (line in ps_output) {
        # Look for typical worker process patterns
        # Workers often have: snow-worker, parallel, SOCK, --max-ppsize, --slave, --no-readline
        if (grepl("snow-worker|parallel|SOCK|--max-ppsize|--slave.*--no-save", line, ignore.case = TRUE)) {
          # Extract PID (2nd column in ps aux output)
          pid <- as.integer(strsplit(trimws(line), "\\s+")[[1]][2])
          
          if (!is.na(pid) && pid != current_pid) {
            kill_cmd <- paste0("kill -9 ", pid, " 2>/dev/null")
            result <- system(kill_cmd)
            
            if (result == 0) {
              killed_processes <- killed_processes + 1
              if (verbose) cat("   ✓ Killed worker process PID:", pid, "\n")
            } else {
              if (verbose) cat("   ✗ Failed to kill PID:", pid, "\n")
            }
          }
        }
      }
      
      if (killed_processes == 0 && verbose) {
        cat("   No worker processes found to kill\n")
      }
    } else {
      if (verbose) cat("   No additional R processes found\n")
    }
  } else {
    if (verbose) {
      cat("\n2. Skipping system process cleanup (kill_all=FALSE)\n")
      cat("   To kill orphaned worker processes, run with kill_all=TRUE\n")
    }
  }
  
  # ============================================================================
  # Step 3: Force garbage collection and report results
  # ============================================================================
  if (verbose) cat("\n3. Running garbage collection...\n")
  
  gc_result <- gc()
  ram_after <- sum(gc_result[,2]) * 1048576
  ram_freed_bytes <- ram_before - ram_after
  ram_freed_gb <- ram_freed_bytes / (1024^3)
  
  if (verbose) {
    cat("   ✓ Garbage collection complete\n")
    cat("\n=== Cleanup Summary ===\n")
    cat("  BiocParallel clusters stopped:", stopped_clusters, "\n")
    cat("  System processes killed:", killed_processes, "\n")
    cat("  RAM freed:", format(ram_freed_bytes, units = "auto", digits = 2), "\n")
    cat("  Current R session RAM:", format(ram_after, units = "auto", digits = 2), "\n")
    cat("\n✓ Cleanup complete!\n\n")
  }
  
  return(invisible(list(
    stopped_clusters = stopped_clusters,
    killed_processes = killed_processes,
    ram_freed_gb = ram_freed_gb,
    ram_current = format(ram_after, units = "auto", digits = 2)
  )))
}

#' Check for Orphaned Workers Without Killing Them
#'
#' Diagnostic function to check for orphaned parallel workers without taking action.
#'
#' @return Character vector of PIDs for potential orphaned workers
#' @export
check_orphaned_workers <- function() {
  cat("\n=== Checking for Orphaned R Workers ===\n")
  
  current_pid <- Sys.getpid()
  ps_cmd <- paste0("ps aux | grep -E 'R|Rscript' | grep -v grep | grep -v '", current_pid, "'")
  ps_output <- system(ps_cmd, intern = TRUE)
  
  if (length(ps_output) == 0) {
    cat("✓ No additional R processes found\n\n")
    return(character(0))
  }
  
  worker_pids <- character()
  
  cat("Found", length(ps_output), "R processes:\n\n")
  for (line in ps_output) {
    # Look for worker indicators
    # Workers have: snow-worker, parallel, SOCK, --max-ppsize, --slave --no-save
    is_worker <- grepl("snow-worker|parallel|SOCK|--max-ppsize|--slave.*--no-save", line, ignore.case = TRUE)
    
    # Extract PID
    pid <- strsplit(trimws(line), "\\s+")[[1]][2]
    
    # Show process with indicator if it's likely a worker
    if (is_worker) {
      cat("⚠️  [WORKER] PID", pid, ":", line, "\n")
      worker_pids <- c(worker_pids, pid)
    } else {
      cat("   [OTHER]  PID", pid, ":", line, "\n")
    }
  }
  
  if (length(worker_pids) > 0) {
    cat("\n⚠️  Found", length(worker_pids), "potential orphaned worker process(es)\n")
    cat("To clean up, run: cleanup_parallel_workers(kill_all = TRUE)\n\n")
  } else {
    cat("\n✓ No orphaned workers detected\n\n")
  }
  
  return(worker_pids)
}

# Auto-run check if sourced directly
if (sys.nframe() == 0) {
  cat("cleanup_parallel_workers.R loaded\n")
  cat("\nUsage:\n")
  cat("  # Safe cleanup (stops registered clusters):\n")
  cat("  cleanup_parallel_workers()\n\n")
  cat("  # Aggressive cleanup (kills all worker processes):\n")
  cat("  cleanup_parallel_workers(kill_all = TRUE)\n\n")
  cat("  # Check for orphaned workers without killing:\n")
  cat("  check_orphaned_workers()\n\n")
}


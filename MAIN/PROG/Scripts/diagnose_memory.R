#' Comprehensive System Memory Diagnostic
#' Author: Dr. Michael Zaiken
#'
#' Checks for memory leaks at multiple levels:
#' - R workspace objects
#' - BiocParallel workers
#' - Java heap (RBioFormats)
#' - Shared memory segments
#' - System process memory
#'
#' @export
diagnose_memory <- function() {
  cat("\n")
  cat("================================================================================\n")
  cat("                    COMPREHENSIVE MEMORY DIAGNOSTIC                             \n")
  cat("================================================================================\n")
  
  # ============================================================================
  # 1. System-Level Memory (Actual RAM Usage)
  # ============================================================================
  cat("\n[1] SYSTEM MEMORY (free -h)\n")
  cat("--------------------------------------------------------------------------------\n")
  system("free -h")
  
  # Get numeric values
  mem_info <- system("free -m | grep '^Mem:'", intern = TRUE)
  mem_parts <- as.numeric(strsplit(trimws(mem_info), "\\s+")[[1]])
  total_ram_mb <- mem_parts[2]
  used_ram_mb <- mem_parts[3]
  
  cat("\n   Total RAM:", total_ram_mb, "MB (", round(total_ram_mb/1024, 1), "GB )\n")
  cat("   Used RAM: ", used_ram_mb, "MB (", round(used_ram_mb/1024, 1), "GB )\n")
  cat("   Usage:    ", round(100 * used_ram_mb / total_ram_mb, 1), "%\n")
  
  # ============================================================================
  # 2. R Session Memory (What R Thinks It's Using)
  # ============================================================================
  cat("\n[2] R SESSION MEMORY (gc())\n")
  cat("--------------------------------------------------------------------------------\n")
  gc_result <- gc()
  print(gc_result)
  r_memory_mb <- sum(gc_result[,2])
  cat("\n   R reports:", round(r_memory_mb, 1), "MB (", round(r_memory_mb/1024, 2), "GB )\n")
  
  # Calculate discrepancy
  discrepancy_mb <- used_ram_mb - r_memory_mb
  if (discrepancy_mb > 1000) {
    cat("   ⚠️  DISCREPANCY:", round(discrepancy_mb, 1), "MB (", round(discrepancy_mb/1024, 1), 
        "GB ) not tracked by R!\n")
  }
  
  # ============================================================================
  # 3. R Workspace Objects
  # ============================================================================
  cat("\n[3] R WORKSPACE OBJECTS (Top 20 by size)\n")
  cat("--------------------------------------------------------------------------------\n")
  obj_sizes <- sapply(ls(envir = .GlobalEnv), function(x) {
    object.size(get(x, envir = .GlobalEnv))
  })
  
  if (length(obj_sizes) > 0) {
    obj_sizes_sorted <- sort(obj_sizes, decreasing = TRUE)
    n_show <- min(20, length(obj_sizes_sorted))
    
    for (i in seq_len(n_show)) {
      cat(sprintf("   %-35s %12s\n", 
                  names(obj_sizes_sorted)[i], 
                  format(obj_sizes_sorted[i], units = "auto")))
    }
    
    total_objects_mb <- sum(obj_sizes) / (1024^2)
    cat("\n   Total workspace objects:", round(total_objects_mb, 1), "MB\n")
  } else {
    cat("   No objects in workspace\n")
  }
  
  # ============================================================================
  # 4. Java Heap Memory (RBioFormats / rJava)
  # ============================================================================
  cat("\n[4] JAVA HEAP MEMORY (if rJava loaded)\n")
  cat("--------------------------------------------------------------------------------\n")
  
  if ("rJava" %in% loadedNamespaces()) {
    tryCatch({
      library(rJava, quietly = TRUE)
      
      # Get Java Runtime
      runtime <- .jcall("java/lang/Runtime", "Ljava/lang/Runtime;", "getRuntime")
      
      # Memory in bytes
      max_mem <- .jcall(runtime, "J", "maxMemory")
      total_mem <- .jcall(runtime, "J", "totalMemory")
      free_mem <- .jcall(runtime, "J", "freeMemory")
      used_mem <- total_mem - free_mem
      
      # Convert to MB
      max_mb <- max_mem / (1024^2)
      total_mb <- total_mem / (1024^2)
      used_mb <- used_mem / (1024^2)
      
      cat("   Java heap max:   ", round(max_mb, 1), "MB (", round(max_mb/1024, 2), "GB )\n")
      cat("   Java heap total: ", round(total_mb, 1), "MB (", round(total_mb/1024, 2), "GB )\n")
      cat("   Java heap used:  ", round(used_mb, 1), "MB (", round(used_mb/1024, 2), "GB )\n")
      
      if (used_mb > 1000) {
        cat("   ⚠️  Java is holding", round(used_mb/1024, 1), "GB!\n")
        cat("   To free: rJava::.jcall('java/lang/System', , 'gc')\n")
      }
    }, error = function(e) {
      cat("   Error checking Java memory:", e$message, "\n")
    })
  } else {
    cat("   rJava not loaded (no Java memory)\n")
  }
  
  # ============================================================================
  # 5. R Process Memory (Actual Process RSS)
  # ============================================================================
  cat("\n[5] R PROCESS MEMORY (ps / top)\n")
  cat("--------------------------------------------------------------------------------\n")
  
  r_pid <- Sys.getpid()
  ps_cmd <- paste0("ps -o pid,rss,vsz,comm -p ", r_pid)
  cat("   PID      RSS(MB)  VSZ(MB)  COMMAND\n")
  ps_output <- system(ps_cmd, intern = TRUE)
  
  if (length(ps_output) > 1) {
    # Parse the output (skip header)
    ps_data <- ps_output[2]
    ps_parts <- strsplit(trimws(ps_data), "\\s+")[[1]]
    
    rss_kb <- as.numeric(ps_parts[2])
    vsz_kb <- as.numeric(ps_parts[3])
    rss_mb <- rss_kb / 1024
    vsz_mb <- vsz_kb / 1024
    
    cat(sprintf("   %-8s %-8.1f %-8.1f %s\n", 
                ps_parts[1], rss_mb, vsz_mb, ps_parts[4]))
    
    cat("\n   R process RSS (Resident Set Size):", round(rss_mb, 1), "MB (", 
        round(rss_mb/1024, 1), "GB )\n")
    
    # This is the actual RAM the R process is using
    if (rss_mb > r_memory_mb * 2) {
      cat("   ⚠️  Process RSS is", round(rss_mb / r_memory_mb, 1), 
          "x larger than R's gc() report!\n")
      cat("   This suggests memory fragmentation or external allocations\n")
    }
  }
  
  # ============================================================================
  # 6. All R Processes (Check for Multiple Sessions / Workers)
  # ============================================================================
  cat("\n[6] ALL R PROCESSES (ps aux | grep R)\n")
  cat("--------------------------------------------------------------------------------\n")
  
  ps_all_cmd <- "ps aux | grep -E '[Rr]session|^.*R |Rscript' | grep -v grep"
  ps_all_output <- system(ps_all_cmd, intern = TRUE)
  
  if (length(ps_all_output) > 0) {
    cat("   Found", length(ps_all_output), "R process(es):\n\n")
    
    total_rss_mb <- 0
    for (line in ps_all_output) {
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      pid <- parts[2]
      rss_pct <- parts[4]  # %MEM
      
      # Get actual RSS for this process
      rss_kb <- as.numeric(system(paste0("ps -p ", pid, " -o rss="), intern = TRUE))
      rss_mb <- rss_kb / 1024
      total_rss_mb <- total_rss_mb + rss_mb
      
      cat(sprintf("   PID %-8s  RSS: %7.1f MB  (%.1f%% MEM)  %s\n", 
                  pid, rss_mb, as.numeric(rss_pct), 
                  substr(paste(parts[11:length(parts)], collapse=" "), 1, 60)))
    }
    
    cat("\n   Total RSS for all R processes:", round(total_rss_mb, 1), 
        "MB (", round(total_rss_mb/1024, 1), "GB )\n")
    
    if (length(ps_all_output) > 1) {
      cat("   ⚠️  Multiple R processes detected! Check for orphaned sessions.\n")
    }
  } else {
    cat("   Only current R process found\n")
  }
  
  # ============================================================================
  # 7. BiocParallel Clusters
  # ============================================================================
  cat("\n[7] BIOCPARALLEL CLUSTERS\n")
  cat("--------------------------------------------------------------------------------\n")
  
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    bpp <- BiocParallel::bpparam()
    cat("   Current backend:", class(bpp)[1], "\n")
    
    if (inherits(bpp, "SnowParam")) {
      n_workers <- BiocParallel::bpnworkers(bpp)
      cat("   ⚠️  ACTIVE SnowParam with", n_workers, "workers detected!\n")
      cat("   Each worker can hold significant memory\n")
      cat("   To stop: parallel::stopCluster(BiocParallel::bpbackend(BiocParallel::bpparam()))\n")
    } else {
      cat("   ✓ No SnowParam cluster (using", class(bpp)[1], ")\n")
    }
  }
  
  # ============================================================================
  # 8. Shared Memory Segments (ipcs)
  # ============================================================================
  cat("\n[8] SHARED MEMORY SEGMENTS (ipcs -m)\n")
  cat("--------------------------------------------------------------------------------\n")
  
  ipcs_output <- system("ipcs -m", intern = TRUE)
  
  # Check if there are any segments
  has_segments <- any(grepl("^0x", ipcs_output))
  
  if (has_segments) {
    cat("   Shared memory segments found:\n\n")
    cat(paste(ipcs_output, collapse = "\n"), "\n")
    
    # Calculate total
    shm_lines <- ipcs_output[grepl("^0x", ipcs_output)]
    if (length(shm_lines) > 0) {
      total_shm_bytes <- sum(sapply(shm_lines, function(line) {
        parts <- strsplit(trimws(line), "\\s+")[[1]]
        as.numeric(parts[5])  # bytes column
      }))
      
      total_shm_mb <- total_shm_bytes / (1024^2)
      if (total_shm_mb > 100) {
        cat("\n   ⚠️  Total shared memory:", round(total_shm_mb, 1), "MB\n")
        cat("   Orphaned shared memory can persist after crashes\n")
      }
    }
  } else {
    cat("   ✓ No shared memory segments\n")
  }
  
  # ============================================================================
  # Summary and Recommendations
  # ============================================================================
  cat("\n")
  cat("================================================================================\n")
  cat("                              SUMMARY                                           \n")
  cat("================================================================================\n")
  
  cat("\n   System RAM used:     ", round(used_ram_mb/1024, 1), "GB\n")
  cat("   R session reports:   ", round(r_memory_mb/1024, 2), "GB\n")
  cat("   Discrepancy:         ", round(discrepancy_mb/1024, 1), "GB\n")
  
  if (discrepancy_mb > 5000) {
    cat("\n   ⚠️⚠️⚠️  LARGE MEMORY LEAK DETECTED (", round(discrepancy_mb/1024, 1), "GB unaccounted) ⚠️⚠️⚠️\n")
    cat("\n   Most likely causes:\n")
    cat("   1. Java heap (RBioFormats) - check section [4]\n")
    cat("   2. Orphaned parallel workers - check section [6]\n")
    cat("   3. Shared memory segments - check section [8]\n")
    cat("   4. Multiple R sessions - check section [6]\n")
    cat("\n   Recommended actions:\n")
    cat("   • Run: cleanup_parallel_workers(kill_all = TRUE)\n")
    cat("   • Run: rJava::.jcall('java/lang/System', , 'gc')\n")
    cat("   • Check for multiple RStudio sessions\n")
    cat("   • Restart R session if issue persists\n")
  }
  
  cat("\n")
  cat("================================================================================\n")
  
  invisible(list(
    system_used_gb = used_ram_mb / 1024,
    r_reports_gb = r_memory_mb / 1024,
    discrepancy_gb = discrepancy_mb / 1024
  ))
}

# Shortcut function to just show the key info
quick_memory_check <- function() {
  cat("\n=== Quick Memory Check ===\n")
  
  # System memory
  mem_info <- system("free -m | grep '^Mem:'", intern = TRUE)
  mem_parts <- as.numeric(strsplit(trimws(mem_info), "\\s+")[[1]])
  system_gb <- mem_parts[3] / 1024
  
  # R memory
  r_gb <- sum(gc()[,2]) / 1024
  
  # R process RSS
  r_pid <- Sys.getpid()
  rss_kb <- as.numeric(system(paste0("ps -p ", r_pid, " -o rss="), intern = TRUE))
  process_gb <- rss_kb / (1024^2)
  
  cat("System used:    ", round(system_gb, 1), "GB\n")
  cat("R session:      ", round(r_gb, 2), "GB\n")
  cat("R process RSS:  ", round(process_gb, 1), "GB\n")
  cat("Leak:           ", round(system_gb - r_gb, 1), "GB\n")
  
  if ((system_gb - r_gb) > 5) {
    cat("\n⚠️  Memory leak detected! Run diagnose_memory() for details\n")
  }
  
  cat("\n")
}

# Auto-display help if sourced
if (sys.nframe() == 0) {
  cat("\nMemory diagnostic tools loaded!\n")
  cat("\nUsage:\n")
  cat("  diagnose_memory()      # Full diagnostic report\n")
  cat("  quick_memory_check()   # Quick summary\n\n")
}


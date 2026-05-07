#!/bin/bash
# Restart R with increased protection stack size
# This fixes "protect(): protection stack overflow" errors

echo "========================================="
echo "Restarting R with stack overflow fixes"
echo "========================================="
echo ""
echo "Settings:"
echo "  --max-ppsize=500000   (protection stack size)"
echo "  R_MAX_VSIZE=256Gb     (max vector heap size)"
echo "  expressions=500000    (max nested expressions, from .Rprofile)"
echo ""
echo "Starting R..."
echo ""

export R_MAX_VSIZE=256Gb
exec R --max-ppsize=500000 --no-save --no-restore

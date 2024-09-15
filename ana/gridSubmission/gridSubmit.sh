# Sample command: . gridSubmit.sh false -14 MnvTunev1_RHC_neg14

neutrinoMode=$1
echo $neutrinoMode

pdg=$2

filename=$3
echo $filename

echo "Let's get started!"
outdir=/pnfs/minerva/persistent/users/mmehmood/default_analysis_loc/$filename

if $neutrinoMode; then 
	echo "Neutrino Mode"
        for STAGE in "eventLoop" "efficiency";do
	#for STAGE in "eventLoop" "migration" "efficiency";do
           for PLAYLIST in "minervame1A" "minervame1B" "minervame1C" "minervame1D" "minervame1E" "minervame1F" "minervame1G" "minervame1L" "minervame1M" "minervame1N" "minervame1O" "minervame1P";do
              python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/mmehmood/MAT_AL9/ --filename $filename --pdg $pdg --outdir $outdir
              echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
           done
        done

   # RUNNING OVER RHC PLAYLISTS
else
	echo "Anti Neutrino mode"
        for STAGE in "eventLoop" "efficiency";do
	#for STAGE in "eventLoop" "migration" "efficiency";do
           for PLAYLIST in "minervame6B" "minervame5A" "minervame6A" "minervame6C" "minervame6D" "minervame6E" "minervame6F" "minervame6G" "minervame6H" "minervame6I" "minervame6J";do
              python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/mmehmood/MAT_AL9/ --filename $filename --pdg $pdg --outdir $outdir
              echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
           done
        done

fi

rm -rf eventLoop_minervame*sh
rm -rf migration_minervame*sh
rm -rf efficiency_minervame*sh 

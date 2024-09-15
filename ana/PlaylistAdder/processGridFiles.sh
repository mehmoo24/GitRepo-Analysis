filename=$1
flagSys=$2

if [ -z "$filename" ]
then
   echo "I'm here to help! Inputs required:"
   echo "Filename (ex. O7cutsO1wgts_nuB_nus)"
   echo "FlagSys (sys or nosys)"
else
   path=/pnfs/minerva/persistent/users/mmehmood/default_analysis_loc/$filename
   cd $path
   for STAGE in "EventSelection" "Migration" "Efficiency";do
      for PLAYLIST in "minervame1A" "minervame1B" "minervame1C" "minervame1D" "minervame1E" "minervame1F" "minervame1G" "minervame1L" "minervame1M" "minervame1N" "minervame1O" "minervame1P";do
         mv Hists_${STAGE}_${filename}_${PLAYLIST}_${flagSys}_t99_z99_Nu.root Hists_${STAGE}_${filename}_${flagSys}_t99_z99_Nu_${PLAYLIST}.root
      done
   done

   cd /exp/minerva/app/users/mmehmood/MAT_AL9/GitRepo-Analysis/ana/PlaylistAdder

   ./PlaylistAdderDebugVersion $path/Hists_EventSelection_${filename}_${flagSys}_t99_z99_Nu_minervame MuonEventSelection 0 1A 1B 1C 1D 1E 1F 1G 1L 1M 1N 1O 1P

   ./PlaylistAdderDebugVersion $path/Hists_Migration_${filename}_${flagSys}_t99_z99_Nu_minervame Migration 0 1A 1B 1C 1D 1E 1F 1G 1L 1M 1N 1O 1P

   ./PlaylistAdderDebugVersion $path/Hists_Efficiency_${filename}_${flagSys}_t99_z99_Nu_minervame EffPurity 0 1A 1B 1C 1D 1E 1F 1G 1L 1M 1N 1O 1P


fi 

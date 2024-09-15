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
      for PLAYLIST in "minervame5A" "minervame6A" "minervame6B" "minervame6C" "minervame6D" "minervame6E" "minervame6F" "minervame6G" "minervame6H" "minervame6I" "minervame6J";do
         mv Hists_${STAGE}_${filename}_${PLAYLIST}_${flagSys}_t99_z99_AntiNu.root Hists_${STAGE}_${filename}_${flagSys}_t99_z99_AntiNu_${PLAYLIST}.root
      done
   done

   cd /exp/minerva/app/users/mmehmood/MAT_AL9/GitRepo-Analysis/ana/PlaylistAdder 
   ./PlaylistAdderDebugVersion $path/Hists_EventSelection_${filename}_${flagSys}_t99_z99_AntiNu_minervame MuonEventSelection 0 5A 6A 6B 6C 6D 6E 6F 6G 6H 6I 6J

   ./PlaylistAdderDebugVersion $path/Hists_Migration_${filename}_${flagSys}_t99_z99_AntiNu_minervame Migration 0 5A 6A 6B 6C 6D 6E 6F 6G 6H 6I 6J

   ./PlaylistAdderDebugVersion $path/Hists_Efficiency_${filename}_${flagSys}_t99_z99_AntiNu_minervame EffPurity 0 5A 6A 6B 6C 6D 6E 6F 6G 6H 6I 6J


fi 

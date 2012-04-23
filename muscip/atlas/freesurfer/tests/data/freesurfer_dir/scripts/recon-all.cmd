#--------------------------------------------
#@# MotionCor Thu Mar 22 18:11:55 EDT 2012

 cp /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/FREESURFER/mri/orig/001.mgz /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/FREESURFER/mri/rawavg.mgz 


 mri_convert /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/FREESURFER/mri/rawavg.mgz /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/FREESURFER/mri/orig.mgz --conform 


 mri_add_xform_to_header -c /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/FREESURFER/mri/transforms/talairach.xfm /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/FREESURFER/mri/orig.mgz /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/FREESURFER/mri/orig.mgz 

#--------------------------------------------
#@# Talairach Thu Mar 22 18:12:08 EDT 2012

 talairach_avi --i orig.mgz --xfm transforms/talairach.auto.xfm 


 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

#--------------------------------------------
#@# Talairach Failure Detection Thu Mar 22 18:12:41 EDT 2012

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /usr/local/freesurfer/bin/extract_talairach_avi_QA.awk /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/FREESURFER/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Nu Intensity Correction Thu Mar 22 18:12:41 EDT 2012

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 

#--------------------------------------------
#@# Intensity Normalization Thu Mar 22 18:14:15 EDT 2012

 mri_normalize -g 1 nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping Thu Mar 22 18:16:39 EDT 2012

 mri_em_register -skull nu.mgz /usr/local/freesurfer/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 


 mri_watershed -T1 -brain_atlas /usr/local/freesurfer/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


 cp brainmask.auto.mgz brainmask.mgz 

#-------------------------------------
#@# EM Registration Thu Mar 22 18:33:20 EDT 2012

 mri_em_register -uns 3 -mask brainmask.mgz nu.mgz /usr/local/freesurfer/average/RB_all_2008-03-26.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize Thu Mar 22 18:51:25 EDT 2012

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /usr/local/freesurfer/average/RB_all_2008-03-26.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg Thu Mar 22 18:53:00 EDT 2012

 mri_ca_register -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /usr/local/freesurfer/average/RB_all_2008-03-26.gca transforms/talairach.m3z 

#--------------------------------------
#@# CA Reg Inv Thu Mar 22 22:39:02 EDT 2012

 mri_ca_register -invert-and-save transforms/talairach.m3z 

#--------------------------------------
#@# Remove Neck Thu Mar 22 22:39:57 EDT 2012

 mri_remove_neck -radius 25 nu.mgz transforms/talairach.m3z /usr/local/freesurfer/average/RB_all_2008-03-26.gca nu_noneck.mgz 

#--------------------------------------
#@# SkullLTA Thu Mar 22 22:40:59 EDT 2012

 mri_em_register -skull -t transforms/talairach.lta nu_noneck.mgz /usr/local/freesurfer/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 

#--------------------------------------
#@# SubCort Seg Thu Mar 22 22:56:50 EDT 2012

 mri_ca_label -align -nobigventricles norm.mgz transforms/talairach.m3z /usr/local/freesurfer/average/RB_all_2008-03-26.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/FREESURFER/mri/transforms/cc_up.lta FREESURFER 

#--------------------------------------
#@# Merge ASeg Thu Mar 22 23:15:35 EDT 2012

 cp aseg.auto.mgz aseg.mgz 

#--------------------------------------------
#@# Intensity Normalization2 Thu Mar 22 23:15:35 EDT 2012

 mri_normalize -aseg aseg.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS Thu Mar 22 23:19:07 EDT 2012

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation Thu Mar 22 23:19:08 EDT 2012

 mri_segment brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.mgz wm.asegedit.mgz 


 mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill Thu Mar 22 23:21:31 EDT 2012

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh Thu Mar 22 23:22:17 EDT 2012

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh Thu Mar 22 23:22:24 EDT 2012

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh Thu Mar 22 23:22:28 EDT 2012

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh Thu Mar 22 23:23:01 EDT 2012

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology lh Thu Mar 22 23:27:18 EDT 2012

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 FREESURFER lh 


 mris_euler_number ../surf/lh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 

#--------------------------------------------
#@# Make White Surf lh Fri Mar 23 00:10:28 EDT 2012

 mris_make_surfaces -noaparc -whiteonly -mgz -T1 brain.finalsurfs FREESURFER lh 

#--------------------------------------------
#@# Smooth2 lh Fri Mar 23 00:15:52 EDT 2012

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white ../surf/lh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh Fri Mar 23 00:15:56 EDT 2012

 mris_inflate ../surf/lh.smoothwm ../surf/lh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/lh.inflated 


#-----------------------------------------
#@# Curvature Stats lh Fri Mar 23 00:18:00 EDT 2012

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm FREESURFER lh curv sulc 

#--------------------------------------------
#@# Sphere lh Fri Mar 23 00:18:05 EDT 2012

 mris_sphere -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Surf Reg lh Fri Mar 23 01:15:43 EDT 2012

 mris_register -curv ../surf/lh.sphere /usr/local/freesurfer/average/lh.average.curvature.filled.buckner40.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh Fri Mar 23 01:46:17 EDT 2012

 mris_jacobian ../surf/lh.white ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh Fri Mar 23 01:46:19 EDT 2012

 mrisp_paint -a 5 /usr/local/freesurfer/average/lh.average.curvature.filled.buckner40.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh Fri Mar 23 01:46:21 EDT 2012

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 FREESURFER lh ../surf/lh.sphere.reg /usr/local/freesurfer/average/lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs ../label/lh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf lh Fri Mar 23 01:47:05 EDT 2012

 mris_make_surfaces -white NOWRITE -mgz -T1 brain.finalsurfs FREESURFER lh 

#--------------------------------------------
#@# Surf Volume lh Fri Mar 23 01:58:17 EDT 2012

 mris_calc -o lh.area.mid lh.area add lh.area.pial 


 mris_calc -o lh.area.mid lh.area.mid div 2 


 mris_calc -o lh.volume lh.area.mid mul lh.thickness 

#-----------------------------------------
#@# Parcellation Stats lh Fri Mar 23 01:58:17 EDT 2012

 mris_anatomical_stats -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab FREESURFER lh white 

#-----------------------------------------
#@# Cortical Parc 2 lh Fri Mar 23 01:58:35 EDT 2012

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 FREESURFER lh ../surf/lh.sphere.reg /usr/local/freesurfer/average/lh.destrieux.simple.2009-07-29.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh Fri Mar 23 01:59:25 EDT 2012

 mris_anatomical_stats -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab FREESURFER lh white 

#--------------------------------------------
#@# Tessellate rh Fri Mar 23 01:59:43 EDT 2012

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 rh Fri Mar 23 01:59:51 EDT 2012

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 rh Fri Mar 23 01:59:55 EDT 2012

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere rh Fri Mar 23 02:00:29 EDT 2012

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology rh Fri Mar 23 02:04:46 EDT 2012

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 FREESURFER rh 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make White Surf rh Fri Mar 23 02:52:18 EDT 2012

 mris_make_surfaces -noaparc -whiteonly -mgz -T1 brain.finalsurfs FREESURFER rh 

#--------------------------------------------
#@# Smooth2 rh Fri Mar 23 02:57:41 EDT 2012

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 rh Fri Mar 23 02:57:45 EDT 2012

 mris_inflate ../surf/rh.smoothwm ../surf/rh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/rh.inflated 


#-----------------------------------------
#@# Curvature Stats rh Fri Mar 23 02:59:48 EDT 2012

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm FREESURFER rh curv sulc 

#--------------------------------------------
#@# Sphere rh Fri Mar 23 02:59:52 EDT 2012

 mris_sphere -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg rh Fri Mar 23 04:01:58 EDT 2012

 mris_register -curv ../surf/rh.sphere /usr/local/freesurfer/average/rh.average.curvature.filled.buckner40.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white rh Fri Mar 23 04:34:01 EDT 2012

 mris_jacobian ../surf/rh.white ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv rh Fri Mar 23 04:34:03 EDT 2012

 mrisp_paint -a 5 /usr/local/freesurfer/average/rh.average.curvature.filled.buckner40.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc rh Fri Mar 23 04:34:05 EDT 2012

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 FREESURFER rh ../surf/rh.sphere.reg /usr/local/freesurfer/average/rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs ../label/rh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf rh Fri Mar 23 04:34:50 EDT 2012

 mris_make_surfaces -white NOWRITE -mgz -T1 brain.finalsurfs FREESURFER rh 

#--------------------------------------------
#@# Surf Volume rh Fri Mar 23 04:45:50 EDT 2012

 mris_calc -o rh.area.mid rh.area add rh.area.pial 


 mris_calc -o rh.area.mid rh.area.mid div 2 


 mris_calc -o rh.volume rh.area.mid mul rh.thickness 

#-----------------------------------------
#@# Parcellation Stats rh Fri Mar 23 04:45:51 EDT 2012

 mris_anatomical_stats -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab FREESURFER rh white 

#-----------------------------------------
#@# Cortical Parc 2 rh Fri Mar 23 04:46:08 EDT 2012

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 FREESURFER rh ../surf/rh.sphere.reg /usr/local/freesurfer/average/rh.destrieux.simple.2009-07-29.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 rh Fri Mar 23 04:46:59 EDT 2012

 mris_anatomical_stats -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab FREESURFER rh white 

#--------------------------------------------
#@# Cortical ribbon mask Fri Mar 23 04:47:18 EDT 2012

 mris_volmask --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon FREESURFER 

#--------------------------------------------
#@# ASeg Stats Fri Mar 23 05:07:50 EDT 2012

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --empty --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --totalgray --ctab /usr/local/freesurfer/ASegStatsLUT.txt --subject FREESURFER 

#-----------------------------------------
#@# AParc-to-ASeg Fri Mar 23 05:15:20 EDT 2012

 mri_aparc2aseg --s FREESURFER --volmask 


 mri_aparc2aseg --s FREESURFER --volmask --a2009s 

#-----------------------------------------
#@# WMParc Fri Mar 23 05:19:38 EDT 2012

 mri_aparc2aseg --s FREESURFER --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brain-vol-from-seg --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --subject FREESURFER --surf-wm-vol --ctab /usr/local/freesurfer/WMParcStatsLUT.txt --etiv 

#--------------------------------------------
#@# BA Labels lh Fri Mar 23 05:34:26 EDT 2012
INFO: fsaverage subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to fsaverage subject...

 cd /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer; ln -s /usr/local/freesurfer/subjects/fsaverage; cd - 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.BA1.label --trgsubject FREESURFER --trglabel ./lh.BA1.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.BA2.label --trgsubject FREESURFER --trglabel ./lh.BA2.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.BA3a.label --trgsubject FREESURFER --trglabel ./lh.BA3a.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.BA3b.label --trgsubject FREESURFER --trglabel ./lh.BA3b.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.BA4a.label --trgsubject FREESURFER --trglabel ./lh.BA4a.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.BA4p.label --trgsubject FREESURFER --trglabel ./lh.BA4p.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.BA6.label --trgsubject FREESURFER --trglabel ./lh.BA6.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.BA44.label --trgsubject FREESURFER --trglabel ./lh.BA44.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.BA45.label --trgsubject FREESURFER --trglabel ./lh.BA45.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.V1.label --trgsubject FREESURFER --trglabel ./lh.V1.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.V2.label --trgsubject FREESURFER --trglabel ./lh.V2.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/lh.MT.label --trgsubject FREESURFER --trglabel ./lh.MT.label --hemi lh --regmethod surface 


 mris_label2annot --s FREESURFER --hemi lh --ctab /usr/local/freesurfer/average/colortable_BA.txt --l lh.BA1.label --l lh.BA2.label --l lh.BA3a.label --l lh.BA3b.label --l lh.BA4a.label --l lh.BA4p.label --l lh.BA6.label --l lh.BA44.label --l lh.BA45.label --l lh.V1.label --l lh.V2.label --l lh.MT.label --a BA --maxstatwinner --noverbose 


 mris_anatomical_stats -mgz -f ../stats/lh.BA.stats -b -a ./lh.BA.annot -c ./BA.ctab FREESURFER lh white 

#--------------------------------------------
#@# BA Labels rh Fri Mar 23 05:36:28 EDT 2012

 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.BA1.label --trgsubject FREESURFER --trglabel ./rh.BA1.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.BA2.label --trgsubject FREESURFER --trglabel ./rh.BA2.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.BA3a.label --trgsubject FREESURFER --trglabel ./rh.BA3a.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.BA3b.label --trgsubject FREESURFER --trglabel ./rh.BA3b.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.BA4a.label --trgsubject FREESURFER --trglabel ./rh.BA4a.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.BA4p.label --trgsubject FREESURFER --trglabel ./rh.BA4p.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.BA6.label --trgsubject FREESURFER --trglabel ./rh.BA6.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.BA44.label --trgsubject FREESURFER --trglabel ./rh.BA44.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.BA45.label --trgsubject FREESURFER --trglabel ./rh.BA45.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.V1.label --trgsubject FREESURFER --trglabel ./rh.V1.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.V2.label --trgsubject FREESURFER --trglabel ./rh.V2.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer/fsaverage/label/rh.MT.label --trgsubject FREESURFER --trglabel ./rh.MT.label --hemi rh --regmethod surface 


 mris_label2annot --s FREESURFER --hemi rh --ctab /usr/local/freesurfer/average/colortable_BA.txt --l rh.BA1.label --l rh.BA2.label --l rh.BA3a.label --l rh.BA3b.label --l rh.BA4a.label --l rh.BA4p.label --l rh.BA6.label --l rh.BA44.label --l rh.BA45.label --l rh.V1.label --l rh.V2.label --l rh.MT.label --a BA --maxstatwinner --noverbose 


 mris_anatomical_stats -mgz -f ../stats/rh.BA.stats -b -a ./rh.BA.annot -c ./BA.ctab FREESURFER rh white 

#--------------------------------------------
#@# Ex-vivo Entorhinal Cortex Label lh Fri Mar 23 05:38:31 EDT 2012
INFO: lh.EC_average subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to lh.EC_average subject...

 cd /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer; ln -s /usr/local/freesurfer/subjects/lh.EC_average; cd - 


 mris_spherical_average -erode 1 -orig white -t 0.4 -o FREESURFER label lh.entorhinal lh sphere.reg lh.EC_average lh.entorhinal_exvivo.label 


 mris_anatomical_stats -mgz -f ../stats/lh.entorhinal_exvivo.stats -b -l ./lh.entorhinal_exvivo.label FREESURFER lh white 

#--------------------------------------------
#@# Ex-vivo Entorhinal Cortex Label rh Fri Mar 23 05:38:44 EDT 2012
INFO: rh.EC_average subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to rh.EC_average subject...

 cd /home/tnesland/Data/tmp/20120321_atlas_compare/NativeFreesurfer; ln -s /usr/local/freesurfer/subjects/rh.EC_average; cd - 


 mris_spherical_average -erode 1 -orig white -t 0.4 -o FREESURFER label rh.entorhinal rh sphere.reg rh.EC_average rh.entorhinal_exvivo.label 


 mris_anatomical_stats -mgz -f ../stats/rh.entorhinal_exvivo.stats -b -l ./rh.entorhinal_exvivo.label FREESURFER rh white 


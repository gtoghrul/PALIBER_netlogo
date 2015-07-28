;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                              ;;
;;                Landscape Fire Succession Model   April 2015                  ;;
;;                                                                              ;;
;;           James D.A. Millington (http://www.landscapemodelling.net)          ;;
;;                                                                              ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;Things to do:
;1) Test (verify) this version of vegetation change
;2) Improve representation of fire (use FireTolerance parameters, frequency, influence on reprouting veg etc)
;3) Use to explore influences of climate vs fire on vegetation change
;4) Add representation of spatially-explicit seed dispersal (currently no seed dispersal represented)
;5) Add representation of spatially-distributed soil conditions (currently spatially uniform across landscape)
;6) Add influence of ungulate (e.g goat) browsing




extensions [ table matrix array ]



globals 
[ 
  pptn-mtx ;matrix containing monthly pptn scenario data (rows = timestep, cols = months)
  temp-mtx ;matrix containing monthly temp scenario data (rows = timestep, cols = months)  
  
  pptn-ary ;array of 12 monthly pptn values
  temp-ary ;array of 12 monthly temp values
 
  PETconstant-ary ;array of constants for 12 months (16 * avegage daylight hrs in month/12 * count days in month/30), does not change through run; from Thornthwaite 1948 eq10 see http://en.wikipedia.org/wiki/Potential_evaporation
  
  
   
  tick-yr ;number of years each tick represents

  fire-front  ;agent-set used in recursive calls in fire simulation
  brow-front  ;agent-set used for browsing simulation
  LCpfs-tab   ;land cover probability of fire spread lookup table used in fire simulation
  BrTol-tab   ;browsing tolerance table used in browsing simulation 
  
  ;filenames 
  filename-LCpfs       ;filename parameter for the land cover probability of fire spread input file
  filename-pptnstream  ;filename parameter for the precipitation weather stream input file
  filename-tempstream  ;filename parameter for the temperature weather stream input file
  filename-pft-mtx     ;filename parameter for the PFT parameter matrix input file
  filename-LCBrTol       ;filename parameter for browsing tolerence (range: 1-3)
  
  pft-mtx  ;matrix holding details of constants for each PFT 
  count-pfts ;the number of pfts being simulated (equal to the length of column 0 of pft-mtx)
  
  DDmin-lst ;list holding Minimum Degree Days for each pft 
  minWT-lst ;list holding mimumum winter temperature for each pft 
  tree-form-lst  ;list holding tree form of each PFT
  maturity-lst   ;list holding maturity ages of each PFT
  resprout-lst   ;list holding resprouting ability of each PFT
  max-age-lst    ;list holding maximum ages of each PFT
  shade-tol-lst    ;list holding shade tolerance of each PFT
  drought-tol-lst  ;list holding drought tolerance of each PFT
  biot-lst         ;list holding form (evergreen or deciduous) of the vegetation
  BioT             ;form of vegetation
  
  kDTT ;development threshold temperature for degree days calc (specified in Table 3.16 of Bugmann 1994 as 5.5)
  elr ;environmental lapse rate (0.00649 degC per 1m, https://en.wikipedia.org/wiki/Lapse_rate#Environmental_lapse_rate)
  
  

]

patches-own 
[
  dlc ;dominant land cover [pft] (the following are always true -2 = unchanging, -1 = bare, then remaining pfts are given by id numbers in pft-mtx)
  prev_lc ;lc in previous tick
  Tin ;time in current state
  prev_Tin ;t_in in previous tick
    
  bucketS ;bucket size  
  PET-pary ;array of 12 monthly PET values, updates each tick
  pptnS-pary ;array of 12 months of pptn surplus (pptn - PET)
  ptemp-pary ;array of 12 months of temperature after environmental lapse rate has been applied
  
  seed-present-plst ;list indicating whether seeds for veg are present
  resp-present-plst ;list indicating whether active resprouter material for veg is present
  established-plst ;list indicating whether a PFT is established in the understory (and therefore can transition to dominant veg) 
  
  age-plst ;age of veg in each PFT (used for tracking maturity of PFT for seed/resprouter material etc.) 
  resources-plst  ;resource level for each PFT 
  
  
  hi ;heat index [can use this to represent influence of temperature on veg? or would that be double accounting as it is used in PET...] 
  alpha ;exponent for PET calc  

  soilM-pary ;array of 12 months of soil moisture 
  Sm-pary ;array of 12 months of S [Bugman & Cramer 1998 Eqn 10] 
  Em-pary ;array of 12 months of E [Bugman & Cramer 1998 Eqn 10] 
  Dm-pary ;array of 12 months of D [Bugman & Cramer 1998 Eqn 10]

  DI ;drought index - calculated from PET 
  uDD ;degree days - calculated from monthly temps
  minT ; minimum winter temperature - calculated from monthly temps
  
  slope 
  aspect  
  SlAsp   ;accounts for slope and aspect in PET (see Bugmann 1994 Eq. 3.73 and 3.74)
  elevation ;elevation of patch (m)
  
  ;slope ;need in futre for fire spread simulation?
  ;aspect ;if slope is reqd for fire spread why not also keep aspect and use some combo to generate SlAsp for each patch in this model?
  
  dD ;direction of transition (land-cover)
  prev_dD ;direction of transition in previous tick
  
  dT ;total time (years) required to complete transition from current to dD
  prev_dT ;time reqd for transition in previous tick
  
  tsf     ;time since fire
  tsb     ;time since browsing
  pfs     ;probability of fire spread
  BrTol   ;browsing tolerance
  BrMor   ;browsing mortality (0 or 1) Bugmann 1994 eq. 3.4, p. 60
  burned  ;did patch burn this time step (use 0/1 for now as may change in future to reflect burn intensity)
  SR ;slope risk used in fire simulation
  FMR ;fuel moisture risk, used in fire simulation
  br-press ;patch specific browsing pressure which decreases with the elevation of the patch
  DrT ;drought tolerance calculated from Henne et al. Table 2
] 


to setup-filenames

  set filename-LCpfs "LC-fire-probs.txt"
  set filename-LCBrTol "LC-BrTol.txt"
  if(temp-scen = "hot") [  set filename-tempstream "wstream-temp-hot.txt" ]
  if(temp-scen = "cool") [  set filename-tempstream "wstream-temp-cool.txt" ]
  if(temp-scen = "average") [ set filename-tempstream "wstream-temp.txt" ]
  
  if(pptn-scen = "wet") [ set filename-pptnstream "wstream-pptn-wet.txt" ]
  if(pptn-scen = "dry") [ set  filename-pptnstream "wstream-pptn-dry.txt" ]
  if(pptn-scen = "average") [ set filename-pptnstream "wstream-pptn.txt" ]
    
  set filename-pft-mtx "pft-parameters_07Apr15.txt"

end

to setup
  
  ca  
  file-close-all
  setup-filenames
  
  reset-ticks
  clear-patches
  
  set kDTT 5.5     ;(specified in Table 3.16 of Bugmann 1994 as 5.5)
  set elr 0.00649  ;environmental lapse rate (degC per m) 
  
  type "setup-relief..."
  setup-relief 
  print "complete"
  
  type "setup-weather-streams..."  
  setup-weather-streams
  print "complete"
  
  type "setup-soil-moisture..."
  setup-soil-moisture
  print "complete"
  
  type "setup-veg..."
  setup-veg
  type "complete (" type count-pfts print " PFTs)"
  
  type "setup-patches..."
  setup-patches
  print "complete"
  
  type "spin-up..."
  ask patches [ spin-up ]  ;execute PET and DI calcs using first year climate data until DI is constant (at given number of dp - currently 5dp)
  print "complete"
  
 
  

  setup-checkParms ;checks parameters max > min etc

  if(fire-sim) 
  [ 
    print "Fire Simulation on"
    setup-SR 
    fire-setup
    setup-FMR
  ]
  
  if (brow-press != 0) [
    print "Browsing Simulation on"
  br-setup
  
  ]
  
  set tick-yr 10 

  clear-all-plots
  reset-ticks
  update-plots
  display-pcolour
      
end

to-report temp-lst
  file-open filename-tempstream
let dummy file-read-line
let this-line []
let out-list []
repeat 12 [
  let temp file-read
  set this-line fput temp this-line
  set out-list fput this-line out-list
  set this-line []
]
file-close
set out-list reverse out-list
report reduce sentence out-list
end

to-report pptn-lst
   file-open filename-pptnstream
let dummy file-read-line
let this-line []
let out-list []
repeat 12 [
  let pptn file-read
  set this-line fput pptn this-line
  set out-list fput this-line out-list
  set this-line []
]
file-close
set out-list reverse out-list
report reduce sentence out-list
end
  


to setup-FMR 

let ann-temp sum temp-lst
let mean-ann-temp ann-temp / 12

let ann-pptn sum pptn-lst
let mean-ann-pptn ann-pptn / 12

let m 12

let sh-par m * (mean-ann-temp / mean-ann-pptn)

ask patches [

ifelse (sh-par < 0.2)
[set FMR 0.8]
[ifelse (sh-par < 0.3)
  [set FMR 0.9]
  [ifelse (sh-par < 0.5)
    [set FMR 1.0]
    [ifelse (sh-par < 0.6)
      [set FMR 1.1]
      [if (sh-par >= 0.6)
        [set FMR 1.2]
      ]
    ]
  ]
]
]
      
end



to setup-SR 
  ask patches [
  ifelse (slope < -25) 
  [set SR 0.80]
  [ifelse (slope < -15) 
    [set SR 0.90]
    [ifelse (slope < -5)
        [set SR 0.95]
        [ifelse (slope < 5)
          [set SR 1.00]
          [ifelse (slope < 15)
            [set SR 1.05]
            [ifelse (slope < 25)
              [set SR 1.10]
              [if (slope > 25)
                [set SR 1.20]
              ]
            ]
          ]
        ]
    ]
  ]
  ]
end



  




to br-setup
  let blist []
  set blist read-BrTol-file filename-LCBrTol blist
  set BrTol-tab table:from-list blist
  type "BrTol:" print BrTol-tab
end



to-report read-BrTol-file [filename outlist]
  file-open filename
  let dummy file-read-line
  
  let this-line []
  let out-list []
  
  while [not file-at-end?]
 [ 
  repeat count-pfts [
    let lc-in file-read
    let pfs-in file-read
    
    set this-line fput pfs-in this-line
    set this-line fput lc-in this-line
    
    set out-list fput this-line out-list
    
    set this-line []  
  ]
 ]
 file-close
 set out-list reverse out-list
 report out-list
end

  





to setup-veg
  
  setup-pft-mtx
  
  set count-pfts length matrix:get-column pft-mtx 0  ;pft-mtx column!
  type "count-pfts... "    ;reports the number of pfts available, 6 are possible at the moment (4 at the spin up) 
  
end
 


to setup-relief

  if(Relief = "Flat") [ ask patches [ set elevation flat-elev ]]
  
  if(Relief = "Valley") [ ask patches [ set elevation abs(pxcor) * 15] ]
  
  if(Relief = "Hilly") 
  [
    ask n-of 30 patches [ set elevation 150000 ] ;just using values that produce feasible elevations (m)
    repeat 1000 [ diffuse elevation 0.525 ]
  ]
  
  if(Relief = "Mountainous") 
  [
    ask n-of 20 patches [ set elevation 800000 ] ;just using values that produce feasible elevations (m)
    repeat 1000 [ diffuse elevation 0.33 ]
  ]    
 
end



  


to setup-checkParms   ;change these to be read from file to matrix holding constants for all PFTs
  
  let column 0
  repeat count-pfts
  [
    if(matrix:get pft-mtx column 1 < 0) [ user-message (sentence "Minimum Degree Days indicator is out of range (<0) in column " column " of " filename-pft-mtx) ] ;pft-mtx column! ;col 1 is minimum degree days
    if(matrix:get pft-mtx column 3 != 0 and matrix:get pft-mtx column 3 != 1) [ user-message (sentence "Tree form indicator is neither 0 nor 1 in column " column " of " filename-pft-mtx) ] ;pft-mtx column! ;col 3 is tree form indicator
    if(matrix:get pft-mtx column 4 < 0 or matrix:get pft-mtx column 4 > 2) [ user-message (sentence "Veg Type (BioT) is out of range in column " column " of " filename-pft-mtx) ] ;pft-mtx column! ;col 4 is BioT
    if(matrix:get pft-mtx column 5 < 0 or matrix:get pft-mtx column 5 > 100) [ user-message (sentence "Maturity age is out of range in column " column " of " filename-pft-mtx) ] ;pft-mtx column! ;col 5 is maturity age
    if(matrix:get pft-mtx column 6 < 0) [ user-message (sentence "Max age is out of range in column " column " of " filename-pft-mtx) ] ;pft-mtx column! ;col 6 is resprout indicator
    if(matrix:get pft-mtx column 7 != 0 and matrix:get pft-mtx column 7 != 1) [ user-message (sentence "Resprout indicator is neither 0 nor 1 in column " column " of " filename-pft-mtx) ] ;pft-mtx column! ;col 7 is resprout indicator
    if(matrix:get pft-mtx column 8 < 1 or matrix:get pft-mtx column 8 > 6) [ user-message (sentence "Shade Tolerance indicator is out of range (<1 or >6) in column " column " of " filename-pft-mtx) ] ;pft-mtx column! ;col 8 is shade tolerance indicator
    if(matrix:get pft-mtx column 9 < 1 or matrix:get pft-mtx column 9 > 6) [ user-message (sentence "Browse Tolerance indicator is out of range (<1 or >6) in column " column " of " filename-pft-mtx) ] ;pft-mtx column! ;col 9 is browse tolerance indicator
    if(matrix:get pft-mtx column 10 < 1 or matrix:get pft-mtx column 10 > 6) [ user-message (sentence "Drought Tolerance indicator is out of range (<1 or >6) in column " column" of " filename-pft-mtx) ] ;pft-mtx column! ;col 10 is drought tolerance indicator
    if(matrix:get pft-mtx column 11 < 1 or matrix:get pft-mtx column 11 > 6) [ user-message (sentence "Fire Tolerance indicator is out of range (<1 or >6) in column " column " of " filename-pft-mtx) ] ;pft-mtx column! ;col 11 is fire tolerance indicator
    
    set column column + 1
  ]
  
  set DDmin-lst matrix:get-column pft-mtx 1 ;pft-mtx column!            ;create list from Minimum Degree Days column of pft matrix
  set minWT-lst matrix:get-column pft-mtx 2 ;pft-mtx column!            ;create list mimumum winter temperature column of pft matrix
  set tree-form-lst matrix:get-column pft-mtx 3  ;pft-mtx column!      ;create list from tree-form column of pft matrix
  set biot-lst matrix:get-column pft-mtx 4     ;create list of BioT column of pft matrix
  set maturity-lst matrix:get-column pft-mtx 5  ;pft-mtx column!      ;create list from maturity column of pft matrix
  set max-age-lst matrix:get-column pft-mtx 6   ;pft-mtx column!      ;create list from maximum age column of pft matrix
  set resprout-lst matrix:get-column pft-mtx 7   ;pft-mtx column!      ;create list from resprouter column of pft matrix
  set shade-tol-lst matrix:get-column pft-mtx 8  ;pft-mtx column!      ;create list from shade tolerance column of pft matrix
  set drought-tol-lst matrix:get-column pft-mtx 10 ;pft-mtx column!      ;create list from drought tolerance column of pft matrix

  
  ;add other parameter checks here
  ;think about how to handle errors when running in batch mode  
   
end


to go

  print date-and-time
 
  ;simulate soil moisture
  ask patches  ;with dlc != -1??
  [
    calc-temp ;calculate temperature for each patch (given environmental lapse rate); include setting minimum winter temp
    calc-PET ;calculate monthly PET
    calc-DI ;calculate drought index
   
    
    calc-uDD ;calculate degree days
    calc-resources ;calculate resources index (from drought and degree days) estab-lst also set here
  ]

  
  ;sim seed dispersal - NOT CURRENTLY IMPLEMENTED 

  
  ;simulate succession  
  update-transitions 
  update-Tin-lc 
  update-vegAge
  
  
  ;sim fire  
  if(fire-sim) [ fire-simulate ]
  if (brow-press != 0) [br-simulate]
  ;output data  
  
  
  display-pcolour
  
  tick
  
end



  



to br-simulate
  if (ticks != 0 and (remainder ticks 5) = 0) [br-start]  ;simulate browsing every 5 steps
end

to br-start
  ask patches with [dlc != -1] [
    set br-press brow-press
    if (elevation > 400) [set br-press brow-press / 2]
    if (elevation > 800) [set br-press brow-press / 4]
    if (elevation > 1200) [set br-press brow-press / 6]
    if (elevation > 1600) [set br-press brow-press / 8]
  set BrTol table:get BrTol-tab dlc ;reads browsing tolerance from the table
  set BrMor (BrTol - 1) * ( (br-press) / 30)  ;Bugmann 1994 eq. 3.4, p.60 
  ]
  if any? patches with [dlc != -1] [ 
  ask patches with [dlc != -1] [     
     if random-float 1 <= BrMor [
       set pcolor blue
       set dlc -1
       set tsb 0
       ]
       ]
       ]
end


   



to fire-simulate
  
  if(ticks != 0 and (remainder ticks 10) = 0) [ fire-ignite ]  ;simulates a fire every 10 steps
  
end

to fire-ignite
  
  ;add code here to determine how frequently fires ignite
  
  ask one-of patches with [dlc != -1]
  [
    set tsf 0
    set dlc -1
    set fire-front patch-set self ;; a new fire-front
    fire-spread
  ]
    
end



to fire-spread

  ;should we continue with the previous method or switch to something more like Henne 2013?
  ;Henne 2013 uses Beta distribution with shape of distro altered (but also acc for drought index) to sim small vs large fires
  ;Henne 2013 approach does not allow representation of influence of individual veg types, slope, wind etc (is that a problem?) but sim is likely much faster 
  

  while [ any? fire-front ]  ;; stop when we run out of fire-front
  [
    show "spreading"
    
    let new-fire-front patch-set nobody ;; empty set of patches for the next 'round' of the fire
    ask fire-front 
    [
      set pcolor red 
      let N neighbors with [ tsf > 0 and dlc != -1]  
      ask N
      [        
        ;calc prob. fire spread here - currently (02Jul13) only uses veg weight (no slope, wind, etc.)
        set pfs table:get LCpfs-tab dlc        
        
        if random-float 1 <= pfs
        [
          set dlc -1
          
          ;what else changes here aside lc? !!CHECK!!
          
          set tsf 0
          set new-fire-front (patch-set new-fire-front self) ;; extend the next round front
        ]
      ]    
    set fire-front new-fire-front
    ]
  ] 
  
end

    

to fire-setup
  
  let flist []
  
  set flist read-LCpfs-file filename-LCpfs flist ;creating list of fire probabilities for different land cover
  
  set LCpfs-tab table:from-list flist ;creating table from the list of fire probabilities for different land cover
  
  type "LC-fire-probs: " print LCpfs-tab
    
end


to-report read-LCpfs-file [filename outlist]
  file-open filename
  let dummy file-read-line ;reads headings to avoid reporting them
  
  let this-line []
  let out-list []
  
  while [not file-at-end?]
  [
    repeat count-pfts 
    [ 
      let lc-in file-read
      let pfs-in file-read
      
      set this-line fput pfs-in this-line
      set this-line fput lc-in this-line
      
      set out-list fput this-line out-list
      
      set this-line[]
    ]
  ]
      
  file-close
  
  set out-list reverse out-list
  
  report out-list
  
end
  

 

to update-vegAge  

;dominant land cover 
;the following are always true -2 = unchanging, -1 = bare, then remaining pfts are given by id numbers)


  ask patches with [burned = 0 and dlc > 0]
  [
    let add-lst[]  ;holds years to add to each pft
    let mature-lst (map [?1 > ?2] age-plst maturity-lst)   ;create list of mature PFTs (if time in state is greater than maturity age)
    
    foreach mature-lst
    [
      ifelse(?) [ set add-lst lput tick-yr add-lst ] [ set add-lst lput 0 add-lst ]  ;convert T/F in mat list to years to add
    ]
    
    set add-lst replace-item dlc add-lst tick-yr  ;overwrite current dlc with tick-yr (don't want to double add for dlc) 
    set age-plst (map + add-lst age-plst)  ;add add-lst to age-plst
  ]
    
    
  ask patches with [burned = 1 and dlc > 0] 
  [
    ;else if burned
    
    let i 0
    foreach resprout-lst 
    [
      ifelse(? = 0) 
      [ set age-plst replace-item i age-plst 0 ]      ;if not a resprouter, kill trees of this PFT (i.e. set age 0)
      [ set seed-present-plst replace-item i seed-present-plst 0 ]    ;if a resprouter, kill seeds of this PFT (i.e. set seed 0) 
      set i i + 1
    ]     
    ;also here check tsf and kill resprouter if necessary? (need a parameter to indicate frequency of burning that kills resprouter material?)
               
  ] 
  
end



;can collapse ask patches code for update-previous into update trans (to save re-asking patches)?
to update-previous
  
  ask patches
  [
    set prev_lc dlc
    set prev_dD dD
    set prev_dT dT
    set prev_Tin Tin
  ]
  
end



;statements 1,2,3,5,6 from Millington et al. (2009)
to update-Tin-lc
  
  ask patches 
  [
    let print-me false
    if(pxcor = 18 and pycor = 15) [ set print-me true ]
    
    if(print-me) [ type "dlc: " type dlc type ", prev_lc: " type prev_lc type ", prev_dD: " type prev_dD type ", prev_dT: " type prev_dT type ", Tin: " print Tin ]
    
    ifelse(dlc != prev_lc) 
    [ set Tin 0 ]
    [
      ifelse(dD != prev_dD) 
      [ set Tin 0 ]                ;lc = prev_lc and dD != prev_dD 
      [ set Tin Tin + tick-yr ]    ;lc = prev_lc and dD = prev_dD
    ]
    
    
    ifelse(Tin > prev_dT) 
    [ 
      set dlc dD
    ]
    [ set dlc prev_lc ]    
    
    
    set tsf tsf + tick-yr   ;code here to update time since fire
    set tsb tsb + tick-yr
         
  ]
  
end
 

to update-transitions
 
  ask patches
  [

    ;for debugging, have only one patch report status
    let print-me false
    if(pxcor = 18 and pycor = 15) [ set print-me true ]

    set prev_lc dlc
    set prev_dD dD
    set prev_dT dT
    set prev_Tin Tin
    
    let trans-set false ;once this is set break out of the procedure
        
    
    ;1. First sort all PFTs on resources (water and growing degree days) descending (H to L)
    let resources-lst-copy resources-plst  ;create a copy of the resources list to identify the order in which PFTs should be checked
    set resources-lst-copy sort-by > resources-lst-copy  ;the sorted resources list (copy)
    
    let positions []  ;this list will tell us the order in which to check the PFTs
    
    foreach resources-lst-copy  ;for each of the sorted list items
    [
      set positions fput position ? resources-plst positions   ;identify the position of the PFT value in the sorted list
    ]
    set positions reverse positions  ; used fput above for speed, so reverse
    ;positions now contains pft id numbers in the order in which they should be checked . 
    
    if(print-me) [ type "positions: " print positions type "resources: " print resources-plst type "established: " print established-plst type "current lc: " print dlc]
    
         
    ;2. Check if highest ranked is current PFT - if so stop... OR this is where we check max age of cohort??
    if(item 0 positions = dlc)
    [
      set dD dlc
      set trans-set true
            
      if(print-me) [ type "dlc: " type dlc type ", position0: " type dD print " set (not checking flags)" ]
    ]
    
    let count-positions length positions
    let this-position 0
    
    ;3. else pick PFT with highest resources in estab-lst (PFT in understory (established-lst true) with max resources)
    if(not trans-set)  ;this will not execute if first PFT had highest resources as trans-set will be true
    [ 
      while[this-position < count-positions]
      [  
        let this-lc item this-position positions
        
        if(item this-lc established-plst = true)
        [
          set dD this-lc
          set trans-set true
          
          if(print-me) [ type "position " type this-position type ": " type dD print " set (not checking flags)" ]
        ]
        
        ifelse(trans-set)
        [ set this-position count-positions if(print-me) [ print "breaking loop"] ] ;will break out of this while loop
        [ set this-position this-position + 1 ] ;will go to next position to check
      ]
    ]  

  
  ;;4 then set time: 10 + dlc-resources * (dD-maturity * (1 + (1 - dD-resources)))
    
    ifelse(dD != dlc)
    [
      let dlc_res 0 ;if dlc < 0 then land is bare so resources are zero 
      if(dlc >= 0) [ set dlc_res item dlc resources-plst ]  ;if dlc >= 0 land is not bare
      
      let u_mat item dD maturity-lst     ;dD is the PFT in the understory transitioning to
      let u_res item dD resources-plst
      
      
      if(print-me) [ type "dlc_res = " type dlc_res type ", u_mat = " type u_mat type ", u_res = " type u_res print "" ]
      
      set dT  precision (10 +  dlc_res * (u_mat * (1 + (1 - u_res)))) -1
      
      if(dT < 0) [ set dT 0 ] ;shouldn't be needed but just in case
      
    ]

    [ set dT 0 ]  ;no time to transition if (dD == dlc)
        
    if(print-me) [ type "dlc = " type dlc type ", dD = " type dD type ", dT = " type dT print "" ]
  ]
  

  
end


  


to calc-resources  ;patch procedure  
  
  ;for debugging, have only one patch report status
  let print-me false
  if(pxcor = 18 and pycor = 15) [ set print-me true ]
     
  set established-plst n-values count-pfts [true]   ;reset estab-lst to all true
 
  
  let column 0
  repeat count-pfts
  [
    let estab true   ;this will change if an establishment flag is set false 
      
    set BioT item column biot-lst
    let DrTolR item column drought-tol-lst
    ifelse (BioT > 0)  [set DrT 0.003 + (DrTolR / 10) * 0.74]
    [set DrT 0.001 + (DrTolR / 10) * 0.82]
    let DDmin item column DDmin-lst
    let minWT item column minWT-lst 
    if(print-me) [ type column type ", DrT: " type DrT type ", DDmin: " type DDmin type ", minT:" type minT type ", minWT: " print minWT ]
    
    let Wres sqrt(max list (1 - (DI / DrT)) 0)   ;water resource; based on Bugmann 1994 Eq. 3.26
     
    let DDres (max list (1 - exp(0.0013 * (DDmin - uDD))) 0) ;degree day resource based on Fyllas et al. 2007 Eq. 9 (and see Fig 2 in Bugmann and Solomon 2000)
        
    if(print-me) [ type column type ", Wres: " type Wres type "DDres: " type DDres ]

        
    ifelse(Wres <= 0 or DDres <= 0)         ;if either resource is unsuitable (i.e. <= 0) set to zero, else calculate average.
    [ set resources-plst replace-item column resources-plst 0 ]     ;position in list corresponds to row in pft-mtx
    [ set resources-plst replace-item column resources-plst ( (Wres + DDres) / 2) ]
     
     
    
     
    ;set estab flag
     ;drought
    if(DI > DrT) [ set estab false ]
    if(print-me) [ type ", DI > DrT: " type DI > DrT ]
               
    ;min temp
    if(minT < minWT) [ set estab false ]
    if(print-me) [ type ", minT < minWT: " type minT < minWT ]
       
    ;degree days
    if(uDD < DDmin) [ set estab false ]
    if(print-me) [ type ", uDD < DDmin: " type uDD < DDmin ]
    
    ;shade tolerance (compare ST of this PFT to ST of other established PFTs)
    let this_ST item column shade-tol-lst
    
    ifelse(dlc < 0) ;if the soil is bare or unchanged, check for shade tolerance.
    [
      if(this_ST > 3) [ set estab false ]
      if(print-me) [ type ", this_ST > 3: " print this_ST > 3 ]  ;assumes PFTs with shade tolerance > 3 cannot establish on bare ground
    ]
    
    [
      let dlc_ST item dlc shade-tol-lst
      if(this_ST < dlc_ST) [ set estab false ] ;comment!!
      if(print-me) [ type ", this_ST < dlc_ST: " print this_ST <= dlc_ST ] ;comparison of shade tolerance of different vegetation
    ]
    
    ;seed/lignotuber availabilty
    
    ;other estab flags?
    
    
    set established-plst replace-item column established-plst estab ;reseting all "true" established-plst with actual "true/false" values
    
    set column column + 1
  ]
  
 
end 



to spin-up 
     

  let last-DI 0  
  let DI-diffc 999
   
  while[DI-diffc != 0]
  [ 
    calc-PET ;calculate monthly PET
    calc-DI ;calculate drought index
    calc-uDD
    
    ;set DI precision DI 3
    set DI-diffc DI - last-DI
    
    set DI-diffc precision DI-diffc 1
    
    set last-DI DI
  ]
  
  ;update-water
  
end


to setup-weather-streams

  let pptn-list []
  set pptn-list read-weather-file filename-pptnstream pptn-list
  set pptn-mtx matrix:from-row-list pptn-list
  ;print matrix:pretty-print-text pptn-mtx

  let temp-list []
  set temp-list read-weather-file filename-tempstream temp-list
  set temp-mtx matrix:from-row-list temp-list
  ;print matrix:pretty-print-text temp-mtx
  
end

to setup-soil-moisture
  
  ;set bucketS soilAWC * soilDepth
  ask patches [
  set bucketS random 21
  if(bucketS < 6) [ set bucketS 6 ]
  if(bucketS > 20) [ set bucketS 20 ]
  set bucketS bucketS * 10 ;converts bucketS from cm to mm (to match PET calcs)
]
  
  let day-hrs [ 9.5  10.75  12  13.5  14.5  15  14.5  13.5  12  10.75  9.5  9.25 ]
  let daylen-ary array:from-list day-hrs

  let mon-days [ 31  28  31  30  31  30  31  31  30  31  30  31 ]
  let monday-ary array:from-list mon-days
  
  let PETconstant []
  let mon 0
  repeat 12
  [
    set PETconstant fput (16 * (array:item daylen-ary mon / 12) * (array:item monday-ary mon / 30) ) PETconstant  ;from Thornthwaite 1948 eq10
    set mon mon + 1
  ]
  
  set PETconstant reverse PETconstant
  set PETconstant-ary array:from-list PETconstant 
  ;print PETconstant-ary   
  
  ask patches
  [
    let empty-list [ ]
    repeat 12 [ set empty-list fput 0 empty-list ]
    set PET-pary array:from-list empty-list
    set pptnS-pary array:from-list empty-list
    
    set soilM-pary array:from-list empty-list 
    ; array:set soilM-ary 11 (matrix:get pptn-mtx 0 0) ;dummy initial value
    
    set Sm-pary array:from-list empty-list
    set Em-pary array:from-list empty-list
    set Dm-pary array:from-list empty-list
  ]
  
    
end



to setup-patches
  
  ask patches
  [
    set dlc -1 ;initially everything is 'bare'
    
    set burned 0
    
    set seed-present-plst n-values count-pfts [true]      ;initially assume all seeds present ;should this list be 0/1 or false/true?
    set resp-present-plst n-values count-pfts [true]      ;initially assume all resprout material present ;should this list be 0/1 or false/true?
    
    set age-plst n-values count-pfts [0]  
    set resources-plst n-values count-pfts [0]  ;should this list be 0/1 or false/true?

    
    ;set dlc random 6
    ;set dlc dlc + 1
    
    set SlAsp 0         ;initially assume SlAsp is always 0 
    
    set dD -1
    set dT tick-yr
    ;update-dDdT
          
    ;set prev_lc -99
    ;set prev_Tin -99
    ;set prev_dD -99
    set prev_dT tick-yr
    
    ;set Tin 0
   
    calc-temp ;calculate temperature for each patch (given environmental lapse rate)

       
  ]

  
end

to-report read-weather-file [ filename wlist ]
 
 file-open filename
  
 let this-line []
 let out-wlist []
 
 let header file-read-line ;read (ignore) the first line which is a header
 
 while [not file-at-end?]
 [
   repeat 12 ;12 months in the year
   [ 
     let in1 file-read
     set this-line fput in1 this-line
   ]
   
   set this-line reverse this-line   ;to report the variables in the correct order
   set out-wlist fput this-line out-wlist
   set this-line []
 ]
 
 file-close
 set out-wlist reverse out-wlist
 report out-wlist
  
end



to setup-pft-mtx
  
  let pft-list []
  set pft-list read-pft-mtx-file filename-pft-mtx pft-list 12 ;currently 12 columns ;pft-mtx column!
  set pft-mtx matrix:from-row-list pft-list
  print ""
  print matrix:pretty-print-text pft-mtx
  
end


to-report read-pft-mtx-file [ filename pft-list columns ]
 
 file-open filename
  
 let this-line []
 let out-pftlist []
 
 let header file-read-line ;read (ignore) the first line which is a header
 
 while [not file-at-end?]
 [
   repeat columns 
   [ 
     let in1 file-read
     set this-line fput in1 this-line
   ]
   
   set this-line reverse this-line
   set out-pftlist fput this-line out-pftlist
   set this-line []
 ]
 
 file-close
 set out-pftlist reverse out-pftlist
 report out-pftlist
  
end




to display-pcolour 
   
  if(Patches-Display = "Vegetation")  ; change this to dominant veg (understory could be different)
  [
    ask patches
    [    
      ifelse(dlc = -1)                 ;bare
      [ set pcolor 4 ]
;      [ ifelse(dlc = 0)               ;grass
;        [ set pcolor 44 ]  
        [ ifelse (dlc = 0)            ;shrubs
          [ set pcolor 24 ]
          [ ifelse (dlc = 1)          ;tree1
            [ set pcolor 74 ]
            [ ifelse (dlc = 2)        ;tree2
              [ set pcolor 54 ]
              [ ifelse (dlc = 3)      ;tree3
                [ set pcolor 34 ]
                [ set pcolor 107 ]  ;water
              ]
            ]
          ]
;        ]
      ]
    ]
  ]
  

 if(Patches-Display = "Drought Index")
 [
   ask patches [ set pcolor scale-color orange DI 0 1 ]
 ]
 
 if(Patches-Display = "Elevation") 
 [ 
   let max-elev 0 
   ask max-one-of patches [elevation] [ set max-elev elevation]
   ask patches [ set pcolor scale-color blue elevation 0 (max-elev + 1) ] ;(max-elev + 1) is used for case if max-elev is 0
 ] 
 
 if(Patches-Display = "Degree Days")
 [
   let max-uDD 0
   ask max-one-of patches [uDD] [ set max-uDD uDD]
   ask patches [ set pcolor scale-color red uDD 0 (max-uDD + 1) ] ; (max-uDD + 1) is used for case if max-uDD is 0
 ] 
 
    
end
      

to calc-temp
  
  let mon 0  ;counter for arrays and matrix columns (i.e. month)
    
  let gtemp-list matrix:get-row temp-mtx ticks ;get monthly landscape temp for this timestep
  let ptemp-list map [? - (elevation * elr)] gtemp-list ;for each month, use elevation of patch with and env lapse rate to modify lsp temp to temp for this patch 
  
  set minT min ptemp-list
    
  ;now edit ptemp-list to not allow negative temps to prevent divide by zero in calc-PET below (hi is sum of ptemp-list and is denominator in eqn in calc-PET). 
  if(min ptemp-list < 0) [ set ptemp-list map [ ifelse-value (? < 0) [0] [?] ] ptemp-list ]  ;if minT is below 0, reports 0; otherwise, reports minT
    
  set ptemp-pary array:from-list ptemp-list  ;set monthly patch temperatures  

end
  

to calc-uDD  ;calculate degree days ;do this once for entire area then modify for individual patches by elevation lapse rate
  
  let mon 0
  
  set uDD 0 
  repeat 12
  [
    let tempt (array:item ptemp-pary mon)
    set tempt tempt - kDTT ;kDTT = 5.5
    if(tempt < 0) [ set tempt 0 ]
    set tempt tempt * 30.5 
    set uDD uDD + tempt 
    set mon mon + 1
  ] 
  
end 


to calc-DI   
  
  ;reminder
  ;soilM-pary ;patch array of 12 months of soil moisture 
  ;Sm-pary ;patch array of 12 months of S [Bugman & Cramer 1998 Eqn 10] 
  ;Em-pary ;patch array of 12 months of E [Bugman & Cramer 1998 Eqn 10] 
  ;Dm-pary ;patch array of 12 months of D [Bugman & Cramer 1998 Eqn 10]
   
  ;calculate soilM for month 0 using final month of last year
  array:set soilM-pary 0 ((array:item soilM-pary 11) + (array:item pptnS-pary 11) - (array:item Em-pary 11))  ; in initial tick all these values will be zero, correct thereafter
  ;Bug & Cram '98 Eq. 9 (soil moisture this month = soil moisture last month + surplus pptn last month - evapotranspiration last month)
  
  if(array:item soilM-pary 0 > bucketS) [ array:set soilM-pary 0 bucketS ] ;if soil moisture is greater than soil capacity, set soil moisture to capcity ;Bug & Cram '98 Eq 9
     ;if soilM > capacity where does the extra water go? runoff downslope??
  

  ;set cw (currently cw is global, but when bucketS becomes patch-based this will be useful here)
  
  let cw 120 ;120mm this is from Henne et al. 2011 Appendix 1 Supplementary material
  if(cw > bucketS) [ set cw bucketS ] 
  
  ;calculate Em for month 0 using last month of last year
  array:set Sm-pary 0 (cw * (array:item soilM-pary 11 / bucketS)) ;Bug & Cram Eq 17 modified by Henne et al. 2011 Appendix 1 Supplementary Material
  array:set Dm-pary 0 (array:item PET-pary 0) ;Bug and Cram '98 Eq. 11 not accounting for intercepted pptn (this will depend on veg type and phenological status)
  
  ifelse(array:item Sm-pary 0 < array:item Dm-pary 0)   ;Bug & Cram 1998 Eq. 10
  [ array:set Em-pary 0 array:item Sm-pary 0 ]  
  [ array:set Em-pary 0 array:item Dm-pary 0 ]
  

  ;now calculate pptnS for all months (pptn supplied in weather stream and PET already calculated for this timestep)
  let mon 0
  repeat 12
  [
    array:set pptnS-pary mon ((matrix:get pptn-mtx ticks mon) - (array:item PET-pary mon))   ;Bugmann and Crammer 1998 Eq. 2
    if(array:item pptnS-pary mon < 0) [ array:set pptnS-pary mon 0 ] 
    set mon mon + 1
  ] 
  
  ;now set soilM-ary for the remaining months (also calculate required Em) 
  set mon 1
  repeat 11
  [
    let last-mon mon - 1
    array:set soilM-pary mon ((array:item soilM-pary last-mon) + (array:item pptnS-pary last-mon) - (array:item Em-pary last-mon))   ;Bug & Cram Eq 9
    if(array:item soilM-pary mon > bucketS) [ array:set soilM-pary mon bucketS ]    ;Bug & Cram Eq 9       ;if soilM > capacity where does the extra water go? runoff downslope??
    
    array:set Sm-pary mon (cw * (array:item soilM-pary last-mon / bucketS))  ;Bug & Cram Eq 17 modified by Henne et al. 2011 Appendix 1 Supplementary Material
    array:set Dm-pary mon (array:item PET-pary mon) ;Bug and Cram '98 Eq. 11 not accounting for intercepted pptn (this will depend on veg type and phenological status)
   
    ifelse(array:item Sm-pary mon < array:item Dm-pary mon)   ;Bug & Cram 1998 Eq. 10
    [ array:set Em-pary mon array:item Sm-pary mon ]  
    [ array:set Em-pary mon array:item Dm-pary mon ]
      
    set mon mon + 1
  ]
  
 
    
  let Em-list array:to-list Em-pary 
  let sum-Em sum Em-list
  
  let Dm-list array:to-list Dm-pary 
  let sum-Dm sum Dm-list
  
  ifelse(sum-Dm = 0) ;may happen because of negative mean monthly temperatures (due to environmental lapse rate) results in PET = 0
  [ set DI 0 ]
  [ set DI 1 - (sum-EM / sum-Dm) ]   ;
  
  if(DI < 0) [ set DI 0 ]   ;do not allow negative DI - !!IS THIS APPROPRIATE/NECESSARY?!!
  
  ;print pptnS-ary
  ;print soilM-ary
  
  ;print Em-ary
  ;print sum-EM
   
  ;print Sm-ary
  ;print sum-PET
  
  ;print DI

end 

 
to calc-PET
  
  ;procedure from Bugmann and Cramer 1996 (drawing on Thornthwaite 1948)
  
  
  let mon 0  ;counter for arrays and matrix columns (i.e. month) 
 
  
  ;calculate heat index and alpha (for PET calculation below)
  let i-list []
  repeat 12 
  [ 
    set i-list fput (((array:item ptemp-pary mon) / 5) ^ 1.514) i-list   ;from Thornthwaite 1948 p89 
    set mon mon + 1    
  ]
  set i-list reverse i-list

  set hi sum i-list
  set alpha (0.000000675 * hi ^ 3) - (0.0000771 * hi ^ 2) + (0.01792 * hi) + 0.49239  ;from Thornthwaite 1948 eq10
  
  ;print i-list
  ;print hi
  ;print alpha
  
  ;all above can be calculated once for entire modelled environment. Only following needs to be done
  
  ;account for slope and aspect see Bugman 1994 Eq 3.73 and 3.74 
  ;let SlAsp 0 ;SlAsp varies by slope and aspect (Steep N slopes = -2, Steep S slopes = +2) - this will need to be incorporated into patch variables
  let kSlAsp 0
  ifelse(SlAsp > 0)
  [ set kSlAsp 1 + (SlAsp * 0.125) ] ;Bugmann 1994 Eq 3.74 - update these parms to better reflect Med conditions?
  [ set kSlAsp 1 + (SlAsp * 0.063) ] ;Bugmann 1994 Eq 3.74 - update these parms to better reflect Med conditions?
    
  
  ;calculate PET by month (note 25May13: confident this is working correctly but not fully tested as currently using same temp etc data for all ticks [i.e. matrix is uniform])
  ;13Jun13 now to accounts for slope and aspect (kSlAsp) 
  
  ifelse ( hi > 0 )  ;hi = 0 [when all months have mean temperature < 0 (as i-list is set to zero in calc-temp procedure)] will cause div/0 error
  [
    set mon 0
    repeat 12
    [
      array:set PET-pary mon ((array:item PETconstant-ary mon) * ((((array:item ptemp-pary mon) * 10) / hi ) ^ alpha )  * kSlAsp )   
      set mon mon + 1
    ]
  ]
  [ ; set PET to 0 for all months
    let zeros n-values 12 [0] 
    set PET-pary array:from-list zeros
  ] 
      
  
end
 

  
@#$#@#$#@
GRAPHICS-WINDOW
210
10
523
344
50
50
3.0
1
10
1
1
1
0
0
0
1
-50
50
-50
50
1
1
1
ticks
30.0

BUTTON
17
22
84
55
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
18
74
81
107
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
17
134
191
167
soilAWC
soilAWC
0.01
0.2
0.2
0.01
1
cm cm-1
HORIZONTAL

SLIDER
18
176
190
209
soilDepth
soilDepth
40
200
80
10
1
cm
HORIZONTAL

BUTTON
92
74
155
107
step
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
412
368
612
518
LC Props
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -5298144 true "" "plot (count patches with [pcolor = 14]) / count patches"
"pen-1" 1.0 0 -15302303 true "" "plot (count patches with [pcolor = 74]) / count patches"
"pen-2" 1.0 0 -14439633 true "" "plot (count patches with [pcolor = 64]) / count patches"
"pen-3" 1.0 0 -12087248 true "" "plot (count patches with [pcolor = 54]) / count patches"
"pen-4" 1.0 0 -8431303 true "" "plot (count patches with [pcolor = 34]) / count patches"
"pen-5" 1.0 0 -3844592 true "" "plot (count patches with [pcolor = 24]) / count patches"

CHOOSER
532
19
673
64
Patches-Display
Patches-Display
"Vegetation" "Drought Index" "Elevation" "Degree Days"
0

BUTTON
536
70
664
103
Update Display
display-pcolour
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
19
220
133
253
fire-sim
fire-sim
0
1
-1000

CHOOSER
18
259
156
304
Relief
Relief
"Flat" "Valley" "Hilly" "Mountainous"
0

PLOT
205
367
405
517
elev-histo
NIL
NIL
0.0
3000.0
0.0
10.0
true
false
"" "histogram [elevation] of patches"
PENS
"pen-0" 1.0 1 -7500403 true "" "histogram [elevation] of patches"

SLIDER
17
313
156
346
flat-elev
flat-elev
0
5000
500
250
1
NIL
HORIZONTAL

CHOOSER
25
374
163
419
temp-scen
temp-scen
"average" "hot" "cool"
0

CHOOSER
24
424
162
469
pptn-scen
pptn-scen
"average" "wet" "dry"
0

SLIDER
647
183
819
216
brow-press
brow-press
0
10
10
0.5
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -7500403 true true 135 285 195 285 270 90 30 90 105 285
Polygon -7500403 true true 270 90 225 15 180 90
Polygon -7500403 true true 30 90 75 15 120 90
Circle -1 true false 183 138 24
Circle -1 true false 93 138 24

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="testing10Dec13" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="40"/>
    <metric>[dlc] of patch 0 0</metric>
    <enumeratedValueSet variable="fire-sim">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="oak-uDD-max">
      <value value="4700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="decid-uDD-min">
      <value value="1300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grass-uDD-max">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pptn-scen">
      <value value="&quot;average&quot;"/>
      <value value="&quot;wet&quot;"/>
      <value value="&quot;dry&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="shrub-uDD-min">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="soilAWC">
      <value value="0.1"/>
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Relief">
      <value value="&quot;Flat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pine-DI-min">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grass-DI-min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flat-elev">
      <value value="0"/>
      <value value="250"/>
      <value value="500"/>
      <value value="750"/>
      <value value="1000"/>
      <value value="1500"/>
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="shrub-uDD-max">
      <value value="6000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="decid-DI-min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pine-uDD-min">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Patches-Display">
      <value value="&quot;Vegetation&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="decid-uDD-max">
      <value value="4700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pine-DI-max">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="oak-uDD-min">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="shrub-DI-max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grass-DI-max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grass-uDD-min">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temp-scen">
      <value value="&quot;average&quot;"/>
      <value value="&quot;hot&quot;"/>
      <value value="&quot;cool&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="decid-DI-max">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pine-uDD-max">
      <value value="2500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Parameters">
      <value value="&quot;From File&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="oak-DI-max">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="oak-DI-min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="shrub-DI-min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="soilDepth">
      <value value="50"/>
      <value value="200"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@

  éI  ½   k820309    ?          18.0        }#ó]                                                                                                          
       src/aed2_phosphorus.F90 AED2_PHOSPHORUS              AED2_PHOSPHORUS_DATA_T                      @                              
                            @                              
       PO4ADSORPTIONFRACTION                   @                               'È                    #AED2_MODEL_ID    #AED2_MODEL_NAME    #AED2_MODEL_PREFIX    #NEXT    #DEFINE    #INITIALIZE    #CALCULATE_SURFACE    #CALCULATE    #CALCULATE_BENTHIC    #CALCULATE_RIPARIAN !   #CALCULATE_DRY '   #EQUILIBRATE ,   #LIGHT_EXTINCTION 1   #RAIN_LOSS 7   #LIGHT_SHADING =   #BIO_DRAG C   #PARTICLE_BGC I   #MOBILITY P   #VALIDATE V   #DELETE [                 $                                                               $                                  @                                   $                                         D                          $                                  È       H             #AED2_MODEL_DATA_T                            ¿              y#AED2_MODEL_DATA_T                                                   1         À    $                                              #AED2_DEFINE 	   #         @     @                            	                    #DATA 
   #NAMLST              
                                
     È               #AED2_MODEL_DATA_T              
                                             1         À    $                                              #AED2_INITIALIZE    #         @     @                                                #DATA    #COLUMN    #LAYER_IDX              
                                      È              #AED2_MODEL_DATA_T              
                                             °                       &                                           #AED2_COLUMN_T              
                                             1         À    $                                              #AED2_CALCULATE_SURFACE    #         @     @                                                #DATA    #COLUMN    #LAYER_IDX              
                                      È              #AED2_MODEL_DATA_T              
                                             °                       &                                           #AED2_COLUMN_T              
                                             1         À    $                                              #AED2_CALCULATE    #         @     @                                                #DATA    #COLUMN    #LAYER_IDX              
                                      È              #AED2_MODEL_DATA_T              
                                             °                       &                                           #AED2_COLUMN_T              
                                             1         À    $                                         	     #AED2_CALCULATE_BENTHIC    #         @     @                                                #DATA    #COLUMN    #LAYER_IDX               
                                      È              #AED2_MODEL_DATA_T              
                                             °                       &                                           #AED2_COLUMN_T              
                                              1         À    $                            !             
     #AED2_CALCULATE_RIPARIAN "   #         @     @                            "                    #DATA #   #COLUMN $   #LAYER_IDX %   #PC_WET &             
                                 #     È              #AED2_MODEL_DATA_T              
                                 $            °       	                &                                           #AED2_COLUMN_T              
                                  %                     
                                 &     
      1         À    $                            '                  #AED2_CALCULATE_DRY (   #         @     @                            (                    #DATA )   #COLUMN *   #LAYER_IDX +             
                                 )     È              #AED2_MODEL_DATA_T              
                                 *            °       
                &                                           #AED2_COLUMN_T              
                                  +           1         À    $                            ,                  #AED2_EQUILIBRATE -   #         @     @                            -                    #DATA .   #COLUMN /   #LAYER_IDX 0             
                                 .     È              #AED2_MODEL_DATA_T              
                                 /            °                       &                                           #AED2_COLUMN_T              
                                  0           1         À    $                            1              	    #AED2_LIGHT_EXTINCTION 2   #         @     @                            2                    #DATA 3   #COLUMN 4   #LAYER_IDX 5   #EXTINCTION 6             
                                 3     È              #AED2_MODEL_DATA_T              
                                 4            °                       &                                           #AED2_COLUMN_T              
                                  5                     
                                6     
       1         À    $                            7              
    #AED2_RAIN_LOSS 8   #         @     @                            8                    #DATA 9   #COLUMN :   #LAYER_IDX ;   #INFIL <             
                                 9     È              #AED2_MODEL_DATA_T              
                                 :            °                       &                                           #AED2_COLUMN_T              
                                  ;                     
                                <     
       1         À    $                            =                  #AED2_LIGHT_SHADING >   #         @     @                            >                    #DATA ?   #COLUMN @   #LAYER_IDX A   #SHADE_FRAC B             
                                 ?     È              #AED2_MODEL_DATA_T              
                                 @            °                       &                                           #AED2_COLUMN_T              
                                  A                     
                                B     
       1         À    $                            C                  #AED2_BIO_DRAG D   #         @     @                            D                    #DATA E   #COLUMN F   #LAYER_IDX G   #DRAG H             
                                 E     È              #AED2_MODEL_DATA_T              
                                 F            °                       &                                           #AED2_COLUMN_T              
                                  G                     
                                H     
       1         À    $                            I                  #AED2_PARTICLE_BGC J   #         @     @                            J                    #DATA K   #COLUMN L   #LAYER_IDX M   #PPID N   #PARTCL O             
                                 K     È              #AED2_MODEL_DATA_T              
                                 L            °                       &                                           #AED2_COLUMN_T              
                                  M                     
                                 N                      
                                O                   
               &                                           1         À    $                            P                  #AED2_MOBILITY Q   #         @     @                            Q                    #DATA R   #COLUMN S   #LAYER_IDX T   #MOBILITY U             
                                 R     È              #AED2_MODEL_DATA_T              
                                 S            °                       &                                           #AED2_COLUMN_T              
                                  T                     
                                U                   
               &                                           1         À    $                           V                  #AED2_VALIDATE W   %         @   @                           W                           #DATA X   #COLUMN Y   #LAYER_IDX Z             
                                 X     È              #AED2_MODEL_DATA_T              
                                 Y            °                       &                                           #AED2_COLUMN_T              
                                  Z           1         À    $                            [                  #AED2_DELETE \   #         @     @                            \                    #DATA ]             
                                ]     È               #AED2_MODEL_DATA_T                      @                                '°                    #CELL ^   #CELL_SHEET _   #FLUX_ATM `   #FLUX_PEL a   #FLUX_BEN b   #FLUX_RIP c               $                             ^                              
            &                                                        $                             _     H          
                $                             `     P          
               $                             a            X                 
            &                                                        $                             b                
                $                             c     ¨          
   #         @                                  d                 	   #PO4ADSORPTIONMODEL e   #PO4TOT_ f   #PARTICLECONC_ g   #KPO4P h   #K i   #QM j   #PO4DIS k   #PO4PAR l   #THEPH m             
                                  e                     
                                 f     
                
                                 g     
                
                                 h     
                
                                 i     
                
                                 j     
                                                k     
                                                 l     
                 
                                m     
                                                  n     
                 
                                 0.                                            o     
                   
                        øÿ                                                    p     
                 
                      õ@        86400.%         @                                q                           #NAME r   #UNITS s   #LONGNAME t   #INITIAL u   #MINIMUM v   #MAXIMUM w   #MOBILITY x             
                                r                    1           
                                s                    1           
                                t                    1           
                                u     
                
                                v     
                
                                w     
                
                                x     
      %         @                                y                           #NAME z             
                                z                    1 %         @                                {                           #NAME |             
                                |                    1 %         @                                }                           #NAME ~             
                                ~                    1 %         @                                                           #NAME    #UNITS    #LONGNAME    #SURF              
                                                    1           
                                                    1           
                                                    1           
                                                              @                              '`             #      #AED2_MODEL_DATA_T    #ID_FRP    #ID_FRPADS    #ID_OXY    #ID_TSS    #ID_PH    #ID_FSED_FRP    #ID_E_TEMP    #ID_E_RAIN    #ID_TSSEXT    #ID_SED_FRP    #ID_FRPADS_VVEL    #ID_ATM_DEP    #FSED_FRP    #KSED_FRP    #THETA_SED_FRP    #ATM_PIP_DD    #ATM_FRP_CONC    #KPO4P    #KADSRATIO    #QMAX    #W_PO4ADS    #SIMDRYDEPOSITION    #SIMWETDEPOSITION    #BEN_USE_OXY    #BEN_USE_AEDSED    #PO4ADSORPTIONMODEL    #SIMPO4ADSORPTION     #ADS_USE_PH ¡   #ADS_USE_EXTERNAL_TSS ¢   #DEFINE £   #CALCULATE_BENTHIC §   #CALCULATE_SURFACE ¬   #EQUILIBRATE ±   #MOBILITY ¶                 $                                   È                      #AED2_MODEL_DATA_T                  $                                   È                           $                                   Ì                           $                                   Ð                           $                                   Ô                           $                                   Ø                           $                                   Ü                           $                                   à                           $                                   ä       	                    $                                   è       
                    $                                   ì                           $                                   ð                           $                                   ô                           $                                  ø          
                 $                                            
                 $                                           
                 $                                           
                 $                                           
                 $                                            
                 $                                  (         
                 $                                  0         
                 $                                  8         
                 $                                   @                          $                                   D                          $                                   H                          $                                   L                          $                                   P                          $                                    T                          $                              ¡     X                          $                              ¢     \            1         À    $                           £                  #AED2_DEFINE_PHOSPHORUS ¤   #         @     @                             ¤                    #DATA ¥   #NAMLST ¦             
D                                ¥     `              #AED2_PHOSPHORUS_DATA_T              
                                  ¦           1         À    $                           §                   #AED2_CALCULATE_BENTHIC_PHOSPHORUS ¨   #         @     @                             ¨                    #DATA ©   #COLUMN ª   #LAYER_IDX «             
                                 ©     `             #AED2_PHOSPHORUS_DATA_T              
D                                 ª            °                       &                                           #AED2_COLUMN_T              
                                  «           1         À    $                           ¬             !     #AED2_CALCULATE_SURFACE_PHOSPHORUS ­   #         @     @                             ­                    #DATA ®   #COLUMN ¯   #LAYER_IDX °             
                                 ®     `             #AED2_PHOSPHORUS_DATA_T              
D                                 ¯            °                       &                                           #AED2_COLUMN_T              
                                  °           1         À    $                           ±             "     #AED2_EQUILIBRATE_PHOSPHORUS ²   #         @     @                             ²                    #DATA ³   #COLUMN ´   #LAYER_IDX µ             
                                 ³     `             #AED2_PHOSPHORUS_DATA_T              
D                                 ´            °                       &                                           #AED2_COLUMN_T              
                                  µ           1         À    $                           ¶             #     #AED2_MOBILITY_PHOSPHORUS ·   #         @     @                             ·                    #DATA ¸   #COLUMN ¹   #LAYER_IDX º   #MOBILITY »             
                                 ¸     `             #AED2_PHOSPHORUS_DATA_T              
                                 ¹            °                       &                                           #AED2_COLUMN_T              
                                  º                     
D                                »                   
               &                                                  0      fn#fn %   Ð   '   b   uapp(AED2_PHOSPHORUS    ÷   @   J  AED2_CORE    7  V   J  AED2_UTIL ,     ®      AED2_MODEL_DATA_T+AED2_CORE :   ;  H   a   AED2_MODEL_DATA_T%AED2_MODEL_ID+AED2_CORE <     P   a   AED2_MODEL_DATA_T%AED2_MODEL_NAME+AED2_CORE >   Ó  P   a   AED2_MODEL_DATA_T%AED2_MODEL_PREFIX+AED2_CORE 1   #  Þ   a   AED2_MODEL_DATA_T%NEXT+AED2_CORE 3     Y   a   AED2_MODEL_DATA_T%DEFINE+AED2_CORE &   Z  ^      AED2_DEFINE+AED2_CORE +   ¸  _   a   AED2_DEFINE%DATA+AED2_CORE -     @   a   AED2_DEFINE%NAMLST+AED2_CORE 7   W  ]   a   AED2_MODEL_DATA_T%INITIALIZE+AED2_CORE *   ´  m      AED2_INITIALIZE+AED2_CORE /   !  _   a   AED2_INITIALIZE%DATA+AED2_CORE 1        a   AED2_INITIALIZE%COLUMN+AED2_CORE 4     @   a   AED2_INITIALIZE%LAYER_IDX+AED2_CORE >   _  d   a   AED2_MODEL_DATA_T%CALCULATE_SURFACE+AED2_CORE 1   Ã  m      AED2_CALCULATE_SURFACE+AED2_CORE 6   0	  _   a   AED2_CALCULATE_SURFACE%DATA+AED2_CORE 8   	     a   AED2_CALCULATE_SURFACE%COLUMN+AED2_CORE ;   .
  @   a   AED2_CALCULATE_SURFACE%LAYER_IDX+AED2_CORE 6   n
  \   a   AED2_MODEL_DATA_T%CALCULATE+AED2_CORE )   Ê
  m      AED2_CALCULATE+AED2_CORE .   7  _   a   AED2_CALCULATE%DATA+AED2_CORE 0        a   AED2_CALCULATE%COLUMN+AED2_CORE 3   5  @   a   AED2_CALCULATE%LAYER_IDX+AED2_CORE >   u  d   a   AED2_MODEL_DATA_T%CALCULATE_BENTHIC+AED2_CORE 1   Ù  m      AED2_CALCULATE_BENTHIC+AED2_CORE 6   F  _   a   AED2_CALCULATE_BENTHIC%DATA+AED2_CORE 8   ¥     a   AED2_CALCULATE_BENTHIC%COLUMN+AED2_CORE ;   D  @   a   AED2_CALCULATE_BENTHIC%LAYER_IDX+AED2_CORE ?     e   a   AED2_MODEL_DATA_T%CALCULATE_RIPARIAN+AED2_CORE 2   é  y      AED2_CALCULATE_RIPARIAN+AED2_CORE 7   b  _   a   AED2_CALCULATE_RIPARIAN%DATA+AED2_CORE 9   Á     a   AED2_CALCULATE_RIPARIAN%COLUMN+AED2_CORE <   `  @   a   AED2_CALCULATE_RIPARIAN%LAYER_IDX+AED2_CORE 9      @   a   AED2_CALCULATE_RIPARIAN%PC_WET+AED2_CORE :   à  `   a   AED2_MODEL_DATA_T%CALCULATE_DRY+AED2_CORE -   @  m      AED2_CALCULATE_DRY+AED2_CORE 2   ­  _   a   AED2_CALCULATE_DRY%DATA+AED2_CORE 4        a   AED2_CALCULATE_DRY%COLUMN+AED2_CORE 7   «  @   a   AED2_CALCULATE_DRY%LAYER_IDX+AED2_CORE 8   ë  ^   a   AED2_MODEL_DATA_T%EQUILIBRATE+AED2_CORE +   I  m      AED2_EQUILIBRATE+AED2_CORE 0   ¶  _   a   AED2_EQUILIBRATE%DATA+AED2_CORE 2        a   AED2_EQUILIBRATE%COLUMN+AED2_CORE 5   ´  @   a   AED2_EQUILIBRATE%LAYER_IDX+AED2_CORE =   ô  c   a   AED2_MODEL_DATA_T%LIGHT_EXTINCTION+AED2_CORE 0   W  }      AED2_LIGHT_EXTINCTION+AED2_CORE 5   Ô  _   a   AED2_LIGHT_EXTINCTION%DATA+AED2_CORE 7   3     a   AED2_LIGHT_EXTINCTION%COLUMN+AED2_CORE :   Ò  @   a   AED2_LIGHT_EXTINCTION%LAYER_IDX+AED2_CORE ;     @   a   AED2_LIGHT_EXTINCTION%EXTINCTION+AED2_CORE 6   R  \   a   AED2_MODEL_DATA_T%RAIN_LOSS+AED2_CORE )   ®  x      AED2_RAIN_LOSS+AED2_CORE .   &  _   a   AED2_RAIN_LOSS%DATA+AED2_CORE 0        a   AED2_RAIN_LOSS%COLUMN+AED2_CORE 3   $  @   a   AED2_RAIN_LOSS%LAYER_IDX+AED2_CORE /   d  @   a   AED2_RAIN_LOSS%INFIL+AED2_CORE :   ¤  `   a   AED2_MODEL_DATA_T%LIGHT_SHADING+AED2_CORE -     }      AED2_LIGHT_SHADING+AED2_CORE 2     _   a   AED2_LIGHT_SHADING%DATA+AED2_CORE 4   à     a   AED2_LIGHT_SHADING%COLUMN+AED2_CORE 7     @   a   AED2_LIGHT_SHADING%LAYER_IDX+AED2_CORE 8   ¿  @   a   AED2_LIGHT_SHADING%SHADE_FRAC+AED2_CORE 5   ÿ  [   a   AED2_MODEL_DATA_T%BIO_DRAG+AED2_CORE (   Z  w      AED2_BIO_DRAG+AED2_CORE -   Ñ  _   a   AED2_BIO_DRAG%DATA+AED2_CORE /   0     a   AED2_BIO_DRAG%COLUMN+AED2_CORE 2   Ï  @   a   AED2_BIO_DRAG%LAYER_IDX+AED2_CORE -     @   a   AED2_BIO_DRAG%DRAG+AED2_CORE 9   O  _   a   AED2_MODEL_DATA_T%PARTICLE_BGC+AED2_CORE ,   ®        AED2_PARTICLE_BGC+AED2_CORE 1   1  _   a   AED2_PARTICLE_BGC%DATA+AED2_CORE 3        a   AED2_PARTICLE_BGC%COLUMN+AED2_CORE 6   /   @   a   AED2_PARTICLE_BGC%LAYER_IDX+AED2_CORE 1   o   @   a   AED2_PARTICLE_BGC%PPID+AED2_CORE 3   ¯      a   AED2_PARTICLE_BGC%PARTCL+AED2_CORE 5   ;!  [   a   AED2_MODEL_DATA_T%MOBILITY+AED2_CORE (   !  {      AED2_MOBILITY+AED2_CORE -   "  _   a   AED2_MOBILITY%DATA+AED2_CORE /   p"     a   AED2_MOBILITY%COLUMN+AED2_CORE 2   #  @   a   AED2_MOBILITY%LAYER_IDX+AED2_CORE 1   O#     a   AED2_MOBILITY%MOBILITY+AED2_CORE 5   Û#  [   a   AED2_MODEL_DATA_T%VALIDATE+AED2_CORE (   6$  u      AED2_VALIDATE+AED2_CORE -   «$  _   a   AED2_VALIDATE%DATA+AED2_CORE /   
%     a   AED2_VALIDATE%COLUMN+AED2_CORE 2   ©%  @   a   AED2_VALIDATE%LAYER_IDX+AED2_CORE 3   é%  Y   a   AED2_MODEL_DATA_T%DELETE+AED2_CORE &   B&  R      AED2_DELETE+AED2_CORE +   &  _   a   AED2_DELETE%DATA+AED2_CORE (   ó&  ¢       AED2_COLUMN_T+AED2_CORE -   '     a   AED2_COLUMN_T%CELL+AED2_CORE 3   )(  H   a   AED2_COLUMN_T%CELL_SHEET+AED2_CORE 1   q(  H   a   AED2_COLUMN_T%FLUX_ATM+AED2_CORE 1   ¹(     a   AED2_COLUMN_T%FLUX_PEL+AED2_CORE 1   M)  H   a   AED2_COLUMN_T%FLUX_BEN+AED2_CORE 1   )  H   a   AED2_COLUMN_T%FLUX_RIP+AED2_CORE 0   Ý)  ½       PO4ADSORPTIONFRACTION+AED2_UTIL C   *  @   a   PO4ADSORPTIONFRACTION%PO4ADSORPTIONMODEL+AED2_UTIL 8   Ú*  @   a   PO4ADSORPTIONFRACTION%PO4TOT_+AED2_UTIL >   +  @   a   PO4ADSORPTIONFRACTION%PARTICLECONC_+AED2_UTIL 6   Z+  @   a   PO4ADSORPTIONFRACTION%KPO4P+AED2_UTIL 2   +  @   a   PO4ADSORPTIONFRACTION%K+AED2_UTIL 3   Ú+  @   a   PO4ADSORPTIONFRACTION%QM+AED2_UTIL 7   ,  @   a   PO4ADSORPTIONFRACTION%PO4DIS+AED2_UTIL 7   Z,  @   a   PO4ADSORPTIONFRACTION%PO4PAR+AED2_UTIL 6   ,  @   a   PO4ADSORPTIONFRACTION%THEPH+AED2_UTIL     Ú,  r       ZERO_+AED2_CORE    L-  p       NAN_+AED2_CORE '   ¼-  v       SECS_PER_DAY+AED2_CORE /   2.  ¨       AED2_DEFINE_VARIABLE+AED2_CORE 4   Ú.  L   a   AED2_DEFINE_VARIABLE%NAME+AED2_CORE 5   &/  L   a   AED2_DEFINE_VARIABLE%UNITS+AED2_CORE 8   r/  L   a   AED2_DEFINE_VARIABLE%LONGNAME+AED2_CORE 7   ¾/  @   a   AED2_DEFINE_VARIABLE%INITIAL+AED2_CORE 7   þ/  @   a   AED2_DEFINE_VARIABLE%MINIMUM+AED2_CORE 7   >0  @   a   AED2_DEFINE_VARIABLE%MAXIMUM+AED2_CORE 8   ~0  @   a   AED2_DEFINE_VARIABLE%MOBILITY+AED2_CORE /   ¾0  Z       AED2_LOCATE_VARIABLE+AED2_CORE 4   1  L   a   AED2_LOCATE_VARIABLE%NAME+AED2_CORE 3   d1  Z       AED2_LOCATE_GLOBAL_SHEET+AED2_CORE 8   ¾1  L   a   AED2_LOCATE_GLOBAL_SHEET%NAME+AED2_CORE -   
2  Z       AED2_LOCATE_GLOBAL+AED2_CORE 2   d2  L   a   AED2_LOCATE_GLOBAL%NAME+AED2_CORE :   °2  }       AED2_DEFINE_SHEET_DIAG_VARIABLE+AED2_CORE ?   -3  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%NAME+AED2_CORE @   y3  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%UNITS+AED2_CORE C   Å3  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%LONGNAME+AED2_CORE ?   4  @   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%SURF+AED2_CORE '   Q4        AED2_PHOSPHORUS_DATA_T 9   í6  g   a   AED2_PHOSPHORUS_DATA_T%AED2_MODEL_DATA_T .   T7  H   a   AED2_PHOSPHORUS_DATA_T%ID_FRP 1   7  H   a   AED2_PHOSPHORUS_DATA_T%ID_FRPADS .   ä7  H   a   AED2_PHOSPHORUS_DATA_T%ID_OXY .   ,8  H   a   AED2_PHOSPHORUS_DATA_T%ID_TSS -   t8  H   a   AED2_PHOSPHORUS_DATA_T%ID_PH 3   ¼8  H   a   AED2_PHOSPHORUS_DATA_T%ID_FSED_FRP 1   9  H   a   AED2_PHOSPHORUS_DATA_T%ID_E_TEMP 1   L9  H   a   AED2_PHOSPHORUS_DATA_T%ID_E_RAIN 1   9  H   a   AED2_PHOSPHORUS_DATA_T%ID_TSSEXT 2   Ü9  H   a   AED2_PHOSPHORUS_DATA_T%ID_SED_FRP 6   $:  H   a   AED2_PHOSPHORUS_DATA_T%ID_FRPADS_VVEL 2   l:  H   a   AED2_PHOSPHORUS_DATA_T%ID_ATM_DEP 0   ´:  H   a   AED2_PHOSPHORUS_DATA_T%FSED_FRP 0   ü:  H   a   AED2_PHOSPHORUS_DATA_T%KSED_FRP 5   D;  H   a   AED2_PHOSPHORUS_DATA_T%THETA_SED_FRP 2   ;  H   a   AED2_PHOSPHORUS_DATA_T%ATM_PIP_DD 4   Ô;  H   a   AED2_PHOSPHORUS_DATA_T%ATM_FRP_CONC -   <  H   a   AED2_PHOSPHORUS_DATA_T%KPO4P 1   d<  H   a   AED2_PHOSPHORUS_DATA_T%KADSRATIO ,   ¬<  H   a   AED2_PHOSPHORUS_DATA_T%QMAX 0   ô<  H   a   AED2_PHOSPHORUS_DATA_T%W_PO4ADS 8   <=  H   a   AED2_PHOSPHORUS_DATA_T%SIMDRYDEPOSITION 8   =  H   a   AED2_PHOSPHORUS_DATA_T%SIMWETDEPOSITION 3   Ì=  H   a   AED2_PHOSPHORUS_DATA_T%BEN_USE_OXY 6   >  H   a   AED2_PHOSPHORUS_DATA_T%BEN_USE_AEDSED :   \>  H   a   AED2_PHOSPHORUS_DATA_T%PO4ADSORPTIONMODEL 8   ¤>  H   a   AED2_PHOSPHORUS_DATA_T%SIMPO4ADSORPTION 2   ì>  H   a   AED2_PHOSPHORUS_DATA_T%ADS_USE_PH <   4?  H   a   AED2_PHOSPHORUS_DATA_T%ADS_USE_EXTERNAL_TSS .   |?  d   a   AED2_PHOSPHORUS_DATA_T%DEFINE '   à?  ^      AED2_DEFINE_PHOSPHORUS ,   >@  d   a   AED2_DEFINE_PHOSPHORUS%DATA .   ¢@  @   a   AED2_DEFINE_PHOSPHORUS%NAMLST 9   â@  o   a   AED2_PHOSPHORUS_DATA_T%CALCULATE_BENTHIC 2   QA  m      AED2_CALCULATE_BENTHIC_PHOSPHORUS 7   ¾A  d   a   AED2_CALCULATE_BENTHIC_PHOSPHORUS%DATA 9   "B     a   AED2_CALCULATE_BENTHIC_PHOSPHORUS%COLUMN <   ÁB  @   a   AED2_CALCULATE_BENTHIC_PHOSPHORUS%LAYER_IDX 9   C  o   a   AED2_PHOSPHORUS_DATA_T%CALCULATE_SURFACE 2   pC  m      AED2_CALCULATE_SURFACE_PHOSPHORUS 7   ÝC  d   a   AED2_CALCULATE_SURFACE_PHOSPHORUS%DATA 9   AD     a   AED2_CALCULATE_SURFACE_PHOSPHORUS%COLUMN <   àD  @   a   AED2_CALCULATE_SURFACE_PHOSPHORUS%LAYER_IDX 3    E  i   a   AED2_PHOSPHORUS_DATA_T%EQUILIBRATE ,   E  m      AED2_EQUILIBRATE_PHOSPHORUS 1   öE  d   a   AED2_EQUILIBRATE_PHOSPHORUS%DATA 3   ZF     a   AED2_EQUILIBRATE_PHOSPHORUS%COLUMN 6   ùF  @   a   AED2_EQUILIBRATE_PHOSPHORUS%LAYER_IDX 0   9G  f   a   AED2_PHOSPHORUS_DATA_T%MOBILITY )   G  {      AED2_MOBILITY_PHOSPHORUS .   H  d   a   AED2_MOBILITY_PHOSPHORUS%DATA 0   ~H     a   AED2_MOBILITY_PHOSPHORUS%COLUMN 3   I  @   a   AED2_MOBILITY_PHOSPHORUS%LAYER_IDX 2   ]I     a   AED2_MOBILITY_PHOSPHORUS%MOBILITY 
  CE  ²   k820309    ?          18.0        |#ó]                                                                                                          
       src/aed2_oxygen.F90 AED2_OXYGEN              AED2_OXYGEN_DATA_T                      @                              
                            @                              
       AED2_GAS_PISTON_VELOCITY AED2_OXYGEN_SAT                   @                               'È                    #AED2_MODEL_ID    #AED2_MODEL_NAME    #AED2_MODEL_PREFIX    #NEXT    #DEFINE    #INITIALIZE    #CALCULATE_SURFACE    #CALCULATE    #CALCULATE_BENTHIC    #CALCULATE_RIPARIAN !   #CALCULATE_DRY '   #EQUILIBRATE ,   #LIGHT_EXTINCTION 1   #RAIN_LOSS 7   #LIGHT_SHADING =   #BIO_DRAG C   #PARTICLE_BGC I   #MOBILITY P   #VALIDATE V   #DELETE [                 $                                                               $                                  @                                   $                                         D                          $                                  È       H             #AED2_MODEL_DATA_T                            9              y#AED2_MODEL_DATA_T                                                   1         À    $                                              #AED2_DEFINE 	   #         @     @                            	                    #DATA 
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
   %         @                              d                 	   
       #WSHGT e   #WIND f   #TEM g   #SAL h   #VEL i   #DEPTH j   #LA k   #SCHMIDT_MODEL l   #PISTON_MODEL m             
                                 e     
                
                                 f     
                
                                 g     
                
                                 h     
                
                                i     
                
                                j     
                
                                k     
                
                                 l                     
                                 m           %         @                              n                    
       #SALT o   #TEMP p             
                                 o     
                
                                 p     
                                                  q     
                   
                        øÿ                                                    r     
                 
                      õ@        86400.%         @                                s                           #NAME t   #UNITS u   #LONGNAME v   #INITIAL w   #MINIMUM x   #MAXIMUM y   #MOBILITY z             
                                t                    1           
                                u                    1           
                                v                    1           
                                w     
                
                                x     
                
                                y     
                
                                z     
      %         @                                {                           #NAME |             
                                |                    1 %         @                                }                           #NAME ~   #UNITS    #LONGNAME              
                                ~                    1           
                                                    1           
                                                    1 %         @                                                           #NAME    #UNITS    #LONGNAME    #SURF              
                                                    1           
                                                    1           
                                                    1           
                                            %         @                                                           #NAME              
                                                    1                                                  
                 
                                 0.                                                 
                 
                       ð?        1.                  @                              '                    #AED2_MODEL_DATA_T    #ID_OXY    #ID_FSED_OXY    #ID_OXY_SAT    #ID_ATM_OXY_EXCH    #ID_SED_OXY    #ID_SED_OXY_PEL    #ID_ATM_OXY_EXCH3D    #ID_TEMP    #ID_SALT    #ID_WIND    #ID_LAREA    #ID_LHT    #ID_CELL_VEL    #FSED_OXY    #KSED_OXY    #THETA_SED_OXY    #OXY_PISTON_MODEL    #USE_SED_MODEL    #DEFINE    #CALCULATE_SURFACE ¢   #CALCULATE §   #CALCULATE_BENTHIC ¬                 $                                   È                      #AED2_MODEL_DATA_T                  $                                   È                           $                                   Ì                           $                                   Ð                           $                                   Ô                           $                                   Ø                           $                                   Ü                           $                                   à                           $                                   ä       	                    $                                   è       
                    $                                   ì                           $                                   ð                           $                                   ô                           $                                   ø                           $                                            
                 $                                           
                 $                                           
                 $                                                             $                                               1         À    $                                             #AED2_DEFINE_OXYGEN    #         @     @                                                 #DATA     #NAMLST ¡             
D                                                     #AED2_OXYGEN_DATA_T              
                                  ¡           1         À    $                           ¢                  #AED2_CALCULATE_SURFACE_OXYGEN £   #         @     @                             £                    #DATA ¤   #COLUMN ¥   #LAYER_IDX ¦             
                                 ¤                   #AED2_OXYGEN_DATA_T              
D                                 ¥            °                       &                                           #AED2_COLUMN_T              
                                  ¦           1         À    $                           §                  #AED2_CALCULATE_OXYGEN ¨   #         @     @                             ¨                    #DATA ©   #COLUMN ª   #LAYER_IDX «             
                                 ©                   #AED2_OXYGEN_DATA_T              
D                                 ª            °                       &                                           #AED2_COLUMN_T              
                                  «           1         À    $                           ¬                  #AED2_CALCULATE_BENTHIC_OXYGEN ­   #         @     @                             ­                    #DATA ®   #COLUMN ¯   #LAYER_IDX °             
                                 ®                   #AED2_OXYGEN_DATA_T              
D                                 ¯            °                       &                                           #AED2_COLUMN_T              
                                  °                  (      fn#fn !   È   #   b   uapp(AED2_OXYGEN    ë   @   J  AED2_CORE    +  i   J  AED2_UTIL ,     ®      AED2_MODEL_DATA_T+AED2_CORE :   B  H   a   AED2_MODEL_DATA_T%AED2_MODEL_ID+AED2_CORE <     P   a   AED2_MODEL_DATA_T%AED2_MODEL_NAME+AED2_CORE >   Ú  P   a   AED2_MODEL_DATA_T%AED2_MODEL_PREFIX+AED2_CORE 1   *  Þ   a   AED2_MODEL_DATA_T%NEXT+AED2_CORE 3     Y   a   AED2_MODEL_DATA_T%DEFINE+AED2_CORE &   a  ^      AED2_DEFINE+AED2_CORE +   ¿  _   a   AED2_DEFINE%DATA+AED2_CORE -     @   a   AED2_DEFINE%NAMLST+AED2_CORE 7   ^  ]   a   AED2_MODEL_DATA_T%INITIALIZE+AED2_CORE *   »  m      AED2_INITIALIZE+AED2_CORE /   (  _   a   AED2_INITIALIZE%DATA+AED2_CORE 1        a   AED2_INITIALIZE%COLUMN+AED2_CORE 4   &  @   a   AED2_INITIALIZE%LAYER_IDX+AED2_CORE >   f  d   a   AED2_MODEL_DATA_T%CALCULATE_SURFACE+AED2_CORE 1   Ê  m      AED2_CALCULATE_SURFACE+AED2_CORE 6   7	  _   a   AED2_CALCULATE_SURFACE%DATA+AED2_CORE 8   	     a   AED2_CALCULATE_SURFACE%COLUMN+AED2_CORE ;   5
  @   a   AED2_CALCULATE_SURFACE%LAYER_IDX+AED2_CORE 6   u
  \   a   AED2_MODEL_DATA_T%CALCULATE+AED2_CORE )   Ñ
  m      AED2_CALCULATE+AED2_CORE .   >  _   a   AED2_CALCULATE%DATA+AED2_CORE 0        a   AED2_CALCULATE%COLUMN+AED2_CORE 3   <  @   a   AED2_CALCULATE%LAYER_IDX+AED2_CORE >   |  d   a   AED2_MODEL_DATA_T%CALCULATE_BENTHIC+AED2_CORE 1   à  m      AED2_CALCULATE_BENTHIC+AED2_CORE 6   M  _   a   AED2_CALCULATE_BENTHIC%DATA+AED2_CORE 8   ¬     a   AED2_CALCULATE_BENTHIC%COLUMN+AED2_CORE ;   K  @   a   AED2_CALCULATE_BENTHIC%LAYER_IDX+AED2_CORE ?     e   a   AED2_MODEL_DATA_T%CALCULATE_RIPARIAN+AED2_CORE 2   ð  y      AED2_CALCULATE_RIPARIAN+AED2_CORE 7   i  _   a   AED2_CALCULATE_RIPARIAN%DATA+AED2_CORE 9   È     a   AED2_CALCULATE_RIPARIAN%COLUMN+AED2_CORE <   g  @   a   AED2_CALCULATE_RIPARIAN%LAYER_IDX+AED2_CORE 9   §  @   a   AED2_CALCULATE_RIPARIAN%PC_WET+AED2_CORE :   ç  `   a   AED2_MODEL_DATA_T%CALCULATE_DRY+AED2_CORE -   G  m      AED2_CALCULATE_DRY+AED2_CORE 2   ´  _   a   AED2_CALCULATE_DRY%DATA+AED2_CORE 4        a   AED2_CALCULATE_DRY%COLUMN+AED2_CORE 7   ²  @   a   AED2_CALCULATE_DRY%LAYER_IDX+AED2_CORE 8   ò  ^   a   AED2_MODEL_DATA_T%EQUILIBRATE+AED2_CORE +   P  m      AED2_EQUILIBRATE+AED2_CORE 0   ½  _   a   AED2_EQUILIBRATE%DATA+AED2_CORE 2        a   AED2_EQUILIBRATE%COLUMN+AED2_CORE 5   »  @   a   AED2_EQUILIBRATE%LAYER_IDX+AED2_CORE =   û  c   a   AED2_MODEL_DATA_T%LIGHT_EXTINCTION+AED2_CORE 0   ^  }      AED2_LIGHT_EXTINCTION+AED2_CORE 5   Û  _   a   AED2_LIGHT_EXTINCTION%DATA+AED2_CORE 7   :     a   AED2_LIGHT_EXTINCTION%COLUMN+AED2_CORE :   Ù  @   a   AED2_LIGHT_EXTINCTION%LAYER_IDX+AED2_CORE ;     @   a   AED2_LIGHT_EXTINCTION%EXTINCTION+AED2_CORE 6   Y  \   a   AED2_MODEL_DATA_T%RAIN_LOSS+AED2_CORE )   µ  x      AED2_RAIN_LOSS+AED2_CORE .   -  _   a   AED2_RAIN_LOSS%DATA+AED2_CORE 0        a   AED2_RAIN_LOSS%COLUMN+AED2_CORE 3   +  @   a   AED2_RAIN_LOSS%LAYER_IDX+AED2_CORE /   k  @   a   AED2_RAIN_LOSS%INFIL+AED2_CORE :   «  `   a   AED2_MODEL_DATA_T%LIGHT_SHADING+AED2_CORE -     }      AED2_LIGHT_SHADING+AED2_CORE 2     _   a   AED2_LIGHT_SHADING%DATA+AED2_CORE 4   ç     a   AED2_LIGHT_SHADING%COLUMN+AED2_CORE 7     @   a   AED2_LIGHT_SHADING%LAYER_IDX+AED2_CORE 8   Æ  @   a   AED2_LIGHT_SHADING%SHADE_FRAC+AED2_CORE 5     [   a   AED2_MODEL_DATA_T%BIO_DRAG+AED2_CORE (   a  w      AED2_BIO_DRAG+AED2_CORE -   Ø  _   a   AED2_BIO_DRAG%DATA+AED2_CORE /   7     a   AED2_BIO_DRAG%COLUMN+AED2_CORE 2   Ö  @   a   AED2_BIO_DRAG%LAYER_IDX+AED2_CORE -     @   a   AED2_BIO_DRAG%DRAG+AED2_CORE 9   V  _   a   AED2_MODEL_DATA_T%PARTICLE_BGC+AED2_CORE ,   µ        AED2_PARTICLE_BGC+AED2_CORE 1   8  _   a   AED2_PARTICLE_BGC%DATA+AED2_CORE 3        a   AED2_PARTICLE_BGC%COLUMN+AED2_CORE 6   6   @   a   AED2_PARTICLE_BGC%LAYER_IDX+AED2_CORE 1   v   @   a   AED2_PARTICLE_BGC%PPID+AED2_CORE 3   ¶      a   AED2_PARTICLE_BGC%PARTCL+AED2_CORE 5   B!  [   a   AED2_MODEL_DATA_T%MOBILITY+AED2_CORE (   !  {      AED2_MOBILITY+AED2_CORE -   "  _   a   AED2_MOBILITY%DATA+AED2_CORE /   w"     a   AED2_MOBILITY%COLUMN+AED2_CORE 2   #  @   a   AED2_MOBILITY%LAYER_IDX+AED2_CORE 1   V#     a   AED2_MOBILITY%MOBILITY+AED2_CORE 5   â#  [   a   AED2_MODEL_DATA_T%VALIDATE+AED2_CORE (   =$  u      AED2_VALIDATE+AED2_CORE -   ²$  _   a   AED2_VALIDATE%DATA+AED2_CORE /   %     a   AED2_VALIDATE%COLUMN+AED2_CORE 2   °%  @   a   AED2_VALIDATE%LAYER_IDX+AED2_CORE 3   ð%  Y   a   AED2_MODEL_DATA_T%DELETE+AED2_CORE &   I&  R      AED2_DELETE+AED2_CORE +   &  _   a   AED2_DELETE%DATA+AED2_CORE (   ú&  ¢       AED2_COLUMN_T+AED2_CORE -   '     a   AED2_COLUMN_T%CELL+AED2_CORE 3   0(  H   a   AED2_COLUMN_T%CELL_SHEET+AED2_CORE 1   x(  H   a   AED2_COLUMN_T%FLUX_ATM+AED2_CORE 1   À(     a   AED2_COLUMN_T%FLUX_PEL+AED2_CORE 1   T)  H   a   AED2_COLUMN_T%FLUX_BEN+AED2_CORE 1   )  H   a   AED2_COLUMN_T%FLUX_RIP+AED2_CORE 3   ä)  ¸       AED2_GAS_PISTON_VELOCITY+AED2_UTIL 9   *  @   a   AED2_GAS_PISTON_VELOCITY%WSHGT+AED2_UTIL 8   Ü*  @   a   AED2_GAS_PISTON_VELOCITY%WIND+AED2_UTIL 7   +  @   a   AED2_GAS_PISTON_VELOCITY%TEM+AED2_UTIL 7   \+  @   a   AED2_GAS_PISTON_VELOCITY%SAL+AED2_UTIL 7   +  @   a   AED2_GAS_PISTON_VELOCITY%VEL+AED2_UTIL 9   Ü+  @   a   AED2_GAS_PISTON_VELOCITY%DEPTH+AED2_UTIL 6   ,  @   a   AED2_GAS_PISTON_VELOCITY%LA+AED2_UTIL A   \,  @   a   AED2_GAS_PISTON_VELOCITY%SCHMIDT_MODEL+AED2_UTIL @   ,  @   a   AED2_GAS_PISTON_VELOCITY%PISTON_MODEL+AED2_UTIL *   Ü,  d       AED2_OXYGEN_SAT+AED2_UTIL /   @-  @   a   AED2_OXYGEN_SAT%SALT+AED2_UTIL /   -  @   a   AED2_OXYGEN_SAT%TEMP+AED2_UTIL    À-  p       NAN_+AED2_CORE '   0.  v       SECS_PER_DAY+AED2_CORE /   ¦.  ¨       AED2_DEFINE_VARIABLE+AED2_CORE 4   N/  L   a   AED2_DEFINE_VARIABLE%NAME+AED2_CORE 5   /  L   a   AED2_DEFINE_VARIABLE%UNITS+AED2_CORE 8   æ/  L   a   AED2_DEFINE_VARIABLE%LONGNAME+AED2_CORE 7   20  @   a   AED2_DEFINE_VARIABLE%INITIAL+AED2_CORE 7   r0  @   a   AED2_DEFINE_VARIABLE%MINIMUM+AED2_CORE 7   ²0  @   a   AED2_DEFINE_VARIABLE%MAXIMUM+AED2_CORE 8   ò0  @   a   AED2_DEFINE_VARIABLE%MOBILITY+AED2_CORE 3   21  Z       AED2_LOCATE_GLOBAL_SHEET+AED2_CORE 8   1  L   a   AED2_LOCATE_GLOBAL_SHEET%NAME+AED2_CORE 4   Ø1  s       AED2_DEFINE_DIAG_VARIABLE+AED2_CORE 9   K2  L   a   AED2_DEFINE_DIAG_VARIABLE%NAME+AED2_CORE :   2  L   a   AED2_DEFINE_DIAG_VARIABLE%UNITS+AED2_CORE =   ã2  L   a   AED2_DEFINE_DIAG_VARIABLE%LONGNAME+AED2_CORE :   /3  }       AED2_DEFINE_SHEET_DIAG_VARIABLE+AED2_CORE ?   ¬3  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%NAME+AED2_CORE @   ø3  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%UNITS+AED2_CORE C   D4  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%LONGNAME+AED2_CORE ?   4  @   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%SURF+AED2_CORE -   Ð4  Z       AED2_LOCATE_GLOBAL+AED2_CORE 2   *5  L   a   AED2_LOCATE_GLOBAL%NAME+AED2_CORE     v5  r       ZERO_+AED2_CORE    è5  r       ONE_+AED2_CORE #   Z6  ×      AED2_OXYGEN_DATA_T 5   18  g   a   AED2_OXYGEN_DATA_T%AED2_MODEL_DATA_T *   8  H   a   AED2_OXYGEN_DATA_T%ID_OXY /   à8  H   a   AED2_OXYGEN_DATA_T%ID_FSED_OXY .   (9  H   a   AED2_OXYGEN_DATA_T%ID_OXY_SAT 3   p9  H   a   AED2_OXYGEN_DATA_T%ID_ATM_OXY_EXCH .   ¸9  H   a   AED2_OXYGEN_DATA_T%ID_SED_OXY 2    :  H   a   AED2_OXYGEN_DATA_T%ID_SED_OXY_PEL 5   H:  H   a   AED2_OXYGEN_DATA_T%ID_ATM_OXY_EXCH3D +   :  H   a   AED2_OXYGEN_DATA_T%ID_TEMP +   Ø:  H   a   AED2_OXYGEN_DATA_T%ID_SALT +    ;  H   a   AED2_OXYGEN_DATA_T%ID_WIND ,   h;  H   a   AED2_OXYGEN_DATA_T%ID_LAREA *   °;  H   a   AED2_OXYGEN_DATA_T%ID_LHT /   ø;  H   a   AED2_OXYGEN_DATA_T%ID_CELL_VEL ,   @<  H   a   AED2_OXYGEN_DATA_T%FSED_OXY ,   <  H   a   AED2_OXYGEN_DATA_T%KSED_OXY 1   Ð<  H   a   AED2_OXYGEN_DATA_T%THETA_SED_OXY 4   =  H   a   AED2_OXYGEN_DATA_T%OXY_PISTON_MODEL 1   `=  H   a   AED2_OXYGEN_DATA_T%USE_SED_MODEL *   ¨=  `   a   AED2_OXYGEN_DATA_T%DEFINE #   >  ^      AED2_DEFINE_OXYGEN (   f>  `   a   AED2_DEFINE_OXYGEN%DATA *   Æ>  @   a   AED2_DEFINE_OXYGEN%NAMLST 5   ?  k   a   AED2_OXYGEN_DATA_T%CALCULATE_SURFACE .   q?  m      AED2_CALCULATE_SURFACE_OXYGEN 3   Þ?  `   a   AED2_CALCULATE_SURFACE_OXYGEN%DATA 5   >@     a   AED2_CALCULATE_SURFACE_OXYGEN%COLUMN 8   Ý@  @   a   AED2_CALCULATE_SURFACE_OXYGEN%LAYER_IDX -   A  c   a   AED2_OXYGEN_DATA_T%CALCULATE &   A  m      AED2_CALCULATE_OXYGEN +   íA  `   a   AED2_CALCULATE_OXYGEN%DATA -   MB     a   AED2_CALCULATE_OXYGEN%COLUMN 0   ìB  @   a   AED2_CALCULATE_OXYGEN%LAYER_IDX 5   ,C  k   a   AED2_OXYGEN_DATA_T%CALCULATE_BENTHIC .   C  m      AED2_CALCULATE_BENTHIC_OXYGEN 3   D  `   a   AED2_CALCULATE_BENTHIC_OXYGEN%DATA 5   dD     a   AED2_CALCULATE_BENTHIC_OXYGEN%COLUMN 8   E  @   a   AED2_CALCULATE_BENTHIC_OXYGEN%LAYER_IDX 
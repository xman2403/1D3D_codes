  -R  ×   k820309    ?          18.0        |#ó]                                                                                                          
       src/aed2_carbon.F90 AED2_CARBON              AED2_CARBON_DATA_T                      @                              
                            @                              
       AED2_GAS_PISTON_VELOCITY                   @                               'È                    #AED2_MODEL_ID    #AED2_MODEL_NAME    #AED2_MODEL_PREFIX    #NEXT    #DEFINE    #INITIALIZE    #CALCULATE_SURFACE    #CALCULATE    #CALCULATE_BENTHIC    #CALCULATE_RIPARIAN !   #CALCULATE_DRY '   #EQUILIBRATE ,   #LIGHT_EXTINCTION 1   #RAIN_LOSS 7   #LIGHT_SHADING =   #BIO_DRAG C   #PARTICLE_BGC I   #MOBILITY P   #VALIDATE V   #DELETE [                 $                                                               $                                  @                                   $                                         D                          $                                  È       H             #AED2_MODEL_DATA_T                                          y#AED2_MODEL_DATA_T                                                   1         À    $                                              #AED2_DEFINE 	   #         @     @                            	                    #DATA 
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
                                 m                                                       n     
                 
                                 0.                                            o     
                 
                       ð?        1.                                            p     
                 
                      õ@        86400.                                            q     
                   
                      ÃÀ        %         @                                r                           #NAME s   #UNITS t   #LONGNAME u   #INITIAL v   #MINIMUM w   #MAXIMUM x   #MOBILITY y             
                                s                    1           
                                t                    1           
                                u                    1           
                                v     
                
                                w     
                
                                x     
                
                                y     
      %         @                                z                           #NAME {             
                                {                    1 %         @                                |                           #NAME }             
                                }                    1 %         @                                ~                           #NAME    #UNITS    #LONGNAME              
                                                    1           
                                                    1           
                                                    1 %         @                                                           #NAME    #UNITS    #LONGNAME    #SURF              
                                                    1           
                                                    1           
                                                    1           
                                            %         @                                                           #NAME              
                                                    1                   @                              'Ø             9      #AED2_MODEL_DATA_T    #ID_DIC    #ID_PH    #ID_CH4    #ID_OXY    #ID_TALK    #ID_CH4_BUB    #ID_FSED_DIC    #ID_FSED_CH4    #ID_TEMP    #ID_SALT    #ID_WIND    #ID_VEL    #ID_DEPTH    #ID_CH4OX    #ID_PCO2    #ID_SED_DIC    #ID_SED_CH4    #ID_SED_CH4_EBB    #ID_ATM_CO2    #ID_ATM_CH4    #ID_ATM_CH4_EBB    #ID_PAR     #ID_EXTC ¡   #ID_DZ ¢   #ID_TAU £   #FSED_DIC ¤   #KSED_DIC ¥   #THETA_SED_DIC ¦   #FSED_CH4 §   #KSED_CH4 ¨   #THETA_SED_CH4 ©   #FSED_CH4_EBB ª   #CH4_BUB_TAU0 «   #RCH4OX ¬   #KCH4OX ­   #VTCH4OX ®   #ATM_CO2 ¯   #ATM_CH4 °   #IONIC ±   #MAXMPBPRODN ²   #IKMPB ³   #USE_OXY ´   #USE_SED_MODEL_DIC µ   #USE_SED_MODEL_CH4 ¶   #SIMDIC ·   #SIMCH4 ¸   #SIMCH4EBB ¹   #ALK_MODE º   #CO2_MODEL »   #CO2_PISTON_MODEL ¼   #CH4_PISTON_MODEL ½   #DEFINE ¾   #CALCULATE_SURFACE Â   #CALCULATE Ç   #CALCULATE_BENTHIC Ì   #EQUILIBRATE Ñ                 $                                   È                      #AED2_MODEL_DATA_T                  $                                   È                           $                                   Ì                           $                                   Ð                           $                                   Ô                           $                                   Ø                           $                                   Ü                           $                                   à                           $                                   ä       	                    $                                   è       
                    $                                   ì                           $                                   ð                           $                                   ô                           $                                   ø                           $                                   ü                           $                                                              $                                                             $                                                             $                                                             $                                                             $                                                             $                                                             $                                                              $                              ¡                                $                              ¢     $                          $                              £     (                          $                             ¤     0         
                 $                             ¥     8         
                 $                             ¦     @         
                 $                             §     H         
                 $                             ¨     P         
                 $                             ©     X          
                 $                             ª     `      !   
                 $                             «     h      "   
                 $                             ¬     p      #   
                 $                             ­     x      $   
                 $                             ®           %   
                 $                             ¯           &   
                 $                             °           '   
                 $                             ±           (   
                 $                             ²            )   
                 $                             ³     ¨      *   
                 $                              ´     °      +                    $                              µ     ´      ,                    $                              ¶     ¸      -                    $                              ·     ¼      .                    $                              ¸     À      /                    $                              ¹     Ä      0                    $                              º     È      1                    $                              »     Ì      2                    $                              ¼     Ð      3                    $                              ½     Ô      4      1         À    $                           ¾             5     #AED2_DEFINE_CARBON ¿   #         @     @                             ¿                    #DATA À   #NAMLST Á             
D                                À     Ø              #AED2_CARBON_DATA_T              
                                  Á           1         À    $                           Â             6     #AED2_CALCULATE_SURFACE_CARBON Ã   #         @     @                             Ã                    #DATA Ä   #COLUMN Å   #LAYER_IDX Æ             
                                 Ä     Ø             #AED2_CARBON_DATA_T              
D                                 Å            °                       &                                           #AED2_COLUMN_T              
                                  Æ           1         À    $                           Ç             7     #AED2_CALCULATE_CARBON È   #         @     @                             È                    #DATA É   #COLUMN Ê   #LAYER_IDX Ë             
                                 É     Ø             #AED2_CARBON_DATA_T              
D                                 Ê            °                       &                                           #AED2_COLUMN_T              
                                  Ë           1         À    $                           Ì             8     #AED2_CALCULATE_BENTHIC_CARBON Í   #         @     @                             Í                    #DATA Î   #COLUMN Ï   #LAYER_IDX Ð             
                                 Î     Ø             #AED2_CARBON_DATA_T              
D                                 Ï            °                       &                                           #AED2_COLUMN_T              
                                  Ð           1         À    $                           Ñ             9     #AED2_EQUILIBRATE_CARBON Ò   #         @     @                             Ò                    #DATA Ó   #COLUMN Ô   #LAYER_IDX Õ             
                                 Ó     Ø             #AED2_CARBON_DATA_T              
D                                 Ô            °                       &                                           #AED2_COLUMN_T              
                                  Õ                  (      fn#fn !   È   #   b   uapp(AED2_CARBON    ë   @   J  AED2_CORE    +  Y   J  AED2_UTIL ,     ®      AED2_MODEL_DATA_T+AED2_CORE :   2  H   a   AED2_MODEL_DATA_T%AED2_MODEL_ID+AED2_CORE <   z  P   a   AED2_MODEL_DATA_T%AED2_MODEL_NAME+AED2_CORE >   Ê  P   a   AED2_MODEL_DATA_T%AED2_MODEL_PREFIX+AED2_CORE 1     Þ   a   AED2_MODEL_DATA_T%NEXT+AED2_CORE 3   ø  Y   a   AED2_MODEL_DATA_T%DEFINE+AED2_CORE &   Q  ^      AED2_DEFINE+AED2_CORE +   ¯  _   a   AED2_DEFINE%DATA+AED2_CORE -     @   a   AED2_DEFINE%NAMLST+AED2_CORE 7   N  ]   a   AED2_MODEL_DATA_T%INITIALIZE+AED2_CORE *   «  m      AED2_INITIALIZE+AED2_CORE /     _   a   AED2_INITIALIZE%DATA+AED2_CORE 1   w     a   AED2_INITIALIZE%COLUMN+AED2_CORE 4     @   a   AED2_INITIALIZE%LAYER_IDX+AED2_CORE >   V  d   a   AED2_MODEL_DATA_T%CALCULATE_SURFACE+AED2_CORE 1   º  m      AED2_CALCULATE_SURFACE+AED2_CORE 6   '	  _   a   AED2_CALCULATE_SURFACE%DATA+AED2_CORE 8   	     a   AED2_CALCULATE_SURFACE%COLUMN+AED2_CORE ;   %
  @   a   AED2_CALCULATE_SURFACE%LAYER_IDX+AED2_CORE 6   e
  \   a   AED2_MODEL_DATA_T%CALCULATE+AED2_CORE )   Á
  m      AED2_CALCULATE+AED2_CORE .   .  _   a   AED2_CALCULATE%DATA+AED2_CORE 0        a   AED2_CALCULATE%COLUMN+AED2_CORE 3   ,  @   a   AED2_CALCULATE%LAYER_IDX+AED2_CORE >   l  d   a   AED2_MODEL_DATA_T%CALCULATE_BENTHIC+AED2_CORE 1   Ð  m      AED2_CALCULATE_BENTHIC+AED2_CORE 6   =  _   a   AED2_CALCULATE_BENTHIC%DATA+AED2_CORE 8        a   AED2_CALCULATE_BENTHIC%COLUMN+AED2_CORE ;   ;  @   a   AED2_CALCULATE_BENTHIC%LAYER_IDX+AED2_CORE ?   {  e   a   AED2_MODEL_DATA_T%CALCULATE_RIPARIAN+AED2_CORE 2   à  y      AED2_CALCULATE_RIPARIAN+AED2_CORE 7   Y  _   a   AED2_CALCULATE_RIPARIAN%DATA+AED2_CORE 9   ¸     a   AED2_CALCULATE_RIPARIAN%COLUMN+AED2_CORE <   W  @   a   AED2_CALCULATE_RIPARIAN%LAYER_IDX+AED2_CORE 9     @   a   AED2_CALCULATE_RIPARIAN%PC_WET+AED2_CORE :   ×  `   a   AED2_MODEL_DATA_T%CALCULATE_DRY+AED2_CORE -   7  m      AED2_CALCULATE_DRY+AED2_CORE 2   ¤  _   a   AED2_CALCULATE_DRY%DATA+AED2_CORE 4        a   AED2_CALCULATE_DRY%COLUMN+AED2_CORE 7   ¢  @   a   AED2_CALCULATE_DRY%LAYER_IDX+AED2_CORE 8   â  ^   a   AED2_MODEL_DATA_T%EQUILIBRATE+AED2_CORE +   @  m      AED2_EQUILIBRATE+AED2_CORE 0   ­  _   a   AED2_EQUILIBRATE%DATA+AED2_CORE 2        a   AED2_EQUILIBRATE%COLUMN+AED2_CORE 5   «  @   a   AED2_EQUILIBRATE%LAYER_IDX+AED2_CORE =   ë  c   a   AED2_MODEL_DATA_T%LIGHT_EXTINCTION+AED2_CORE 0   N  }      AED2_LIGHT_EXTINCTION+AED2_CORE 5   Ë  _   a   AED2_LIGHT_EXTINCTION%DATA+AED2_CORE 7   *     a   AED2_LIGHT_EXTINCTION%COLUMN+AED2_CORE :   É  @   a   AED2_LIGHT_EXTINCTION%LAYER_IDX+AED2_CORE ;   	  @   a   AED2_LIGHT_EXTINCTION%EXTINCTION+AED2_CORE 6   I  \   a   AED2_MODEL_DATA_T%RAIN_LOSS+AED2_CORE )   ¥  x      AED2_RAIN_LOSS+AED2_CORE .     _   a   AED2_RAIN_LOSS%DATA+AED2_CORE 0   |     a   AED2_RAIN_LOSS%COLUMN+AED2_CORE 3     @   a   AED2_RAIN_LOSS%LAYER_IDX+AED2_CORE /   [  @   a   AED2_RAIN_LOSS%INFIL+AED2_CORE :     `   a   AED2_MODEL_DATA_T%LIGHT_SHADING+AED2_CORE -   û  }      AED2_LIGHT_SHADING+AED2_CORE 2   x  _   a   AED2_LIGHT_SHADING%DATA+AED2_CORE 4   ×     a   AED2_LIGHT_SHADING%COLUMN+AED2_CORE 7   v  @   a   AED2_LIGHT_SHADING%LAYER_IDX+AED2_CORE 8   ¶  @   a   AED2_LIGHT_SHADING%SHADE_FRAC+AED2_CORE 5   ö  [   a   AED2_MODEL_DATA_T%BIO_DRAG+AED2_CORE (   Q  w      AED2_BIO_DRAG+AED2_CORE -   È  _   a   AED2_BIO_DRAG%DATA+AED2_CORE /   '     a   AED2_BIO_DRAG%COLUMN+AED2_CORE 2   Æ  @   a   AED2_BIO_DRAG%LAYER_IDX+AED2_CORE -     @   a   AED2_BIO_DRAG%DRAG+AED2_CORE 9   F  _   a   AED2_MODEL_DATA_T%PARTICLE_BGC+AED2_CORE ,   ¥        AED2_PARTICLE_BGC+AED2_CORE 1   (  _   a   AED2_PARTICLE_BGC%DATA+AED2_CORE 3        a   AED2_PARTICLE_BGC%COLUMN+AED2_CORE 6   &   @   a   AED2_PARTICLE_BGC%LAYER_IDX+AED2_CORE 1   f   @   a   AED2_PARTICLE_BGC%PPID+AED2_CORE 3   ¦      a   AED2_PARTICLE_BGC%PARTCL+AED2_CORE 5   2!  [   a   AED2_MODEL_DATA_T%MOBILITY+AED2_CORE (   !  {      AED2_MOBILITY+AED2_CORE -   "  _   a   AED2_MOBILITY%DATA+AED2_CORE /   g"     a   AED2_MOBILITY%COLUMN+AED2_CORE 2   #  @   a   AED2_MOBILITY%LAYER_IDX+AED2_CORE 1   F#     a   AED2_MOBILITY%MOBILITY+AED2_CORE 5   Ò#  [   a   AED2_MODEL_DATA_T%VALIDATE+AED2_CORE (   -$  u      AED2_VALIDATE+AED2_CORE -   ¢$  _   a   AED2_VALIDATE%DATA+AED2_CORE /   %     a   AED2_VALIDATE%COLUMN+AED2_CORE 2    %  @   a   AED2_VALIDATE%LAYER_IDX+AED2_CORE 3   à%  Y   a   AED2_MODEL_DATA_T%DELETE+AED2_CORE &   9&  R      AED2_DELETE+AED2_CORE +   &  _   a   AED2_DELETE%DATA+AED2_CORE (   ê&  ¢       AED2_COLUMN_T+AED2_CORE -   '     a   AED2_COLUMN_T%CELL+AED2_CORE 3    (  H   a   AED2_COLUMN_T%CELL_SHEET+AED2_CORE 1   h(  H   a   AED2_COLUMN_T%FLUX_ATM+AED2_CORE 1   °(     a   AED2_COLUMN_T%FLUX_PEL+AED2_CORE 1   D)  H   a   AED2_COLUMN_T%FLUX_BEN+AED2_CORE 1   )  H   a   AED2_COLUMN_T%FLUX_RIP+AED2_CORE 3   Ô)  ¸       AED2_GAS_PISTON_VELOCITY+AED2_UTIL 9   *  @   a   AED2_GAS_PISTON_VELOCITY%WSHGT+AED2_UTIL 8   Ì*  @   a   AED2_GAS_PISTON_VELOCITY%WIND+AED2_UTIL 7   +  @   a   AED2_GAS_PISTON_VELOCITY%TEM+AED2_UTIL 7   L+  @   a   AED2_GAS_PISTON_VELOCITY%SAL+AED2_UTIL 7   +  @   a   AED2_GAS_PISTON_VELOCITY%VEL+AED2_UTIL 9   Ì+  @   a   AED2_GAS_PISTON_VELOCITY%DEPTH+AED2_UTIL 6   ,  @   a   AED2_GAS_PISTON_VELOCITY%LA+AED2_UTIL A   L,  @   a   AED2_GAS_PISTON_VELOCITY%SCHMIDT_MODEL+AED2_UTIL @   ,  @   a   AED2_GAS_PISTON_VELOCITY%PISTON_MODEL+AED2_UTIL     Ì,  r       ZERO_+AED2_CORE    >-  r       ONE_+AED2_CORE '   °-  v       SECS_PER_DAY+AED2_CORE "   &.  p       MISVAL_+AED2_CORE /   .  ¨       AED2_DEFINE_VARIABLE+AED2_CORE 4   >/  L   a   AED2_DEFINE_VARIABLE%NAME+AED2_CORE 5   /  L   a   AED2_DEFINE_VARIABLE%UNITS+AED2_CORE 8   Ö/  L   a   AED2_DEFINE_VARIABLE%LONGNAME+AED2_CORE 7   "0  @   a   AED2_DEFINE_VARIABLE%INITIAL+AED2_CORE 7   b0  @   a   AED2_DEFINE_VARIABLE%MINIMUM+AED2_CORE 7   ¢0  @   a   AED2_DEFINE_VARIABLE%MAXIMUM+AED2_CORE 8   â0  @   a   AED2_DEFINE_VARIABLE%MOBILITY+AED2_CORE /   "1  Z       AED2_LOCATE_VARIABLE+AED2_CORE 4   |1  L   a   AED2_LOCATE_VARIABLE%NAME+AED2_CORE 3   È1  Z       AED2_LOCATE_GLOBAL_SHEET+AED2_CORE 8   "2  L   a   AED2_LOCATE_GLOBAL_SHEET%NAME+AED2_CORE 4   n2  s       AED2_DEFINE_DIAG_VARIABLE+AED2_CORE 9   á2  L   a   AED2_DEFINE_DIAG_VARIABLE%NAME+AED2_CORE :   -3  L   a   AED2_DEFINE_DIAG_VARIABLE%UNITS+AED2_CORE =   y3  L   a   AED2_DEFINE_DIAG_VARIABLE%LONGNAME+AED2_CORE :   Å3  }       AED2_DEFINE_SHEET_DIAG_VARIABLE+AED2_CORE ?   B4  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%NAME+AED2_CORE @   4  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%UNITS+AED2_CORE C   Ú4  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%LONGNAME+AED2_CORE ?   &5  @   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%SURF+AED2_CORE -   f5  Z       AED2_LOCATE_GLOBAL+AED2_CORE 2   À5  L   a   AED2_LOCATE_GLOBAL%NAME+AED2_CORE #   6  ¶      AED2_CARBON_DATA_T 5   Â9  g   a   AED2_CARBON_DATA_T%AED2_MODEL_DATA_T *   ):  H   a   AED2_CARBON_DATA_T%ID_DIC )   q:  H   a   AED2_CARBON_DATA_T%ID_PH *   ¹:  H   a   AED2_CARBON_DATA_T%ID_CH4 *   ;  H   a   AED2_CARBON_DATA_T%ID_OXY +   I;  H   a   AED2_CARBON_DATA_T%ID_TALK .   ;  H   a   AED2_CARBON_DATA_T%ID_CH4_BUB /   Ù;  H   a   AED2_CARBON_DATA_T%ID_FSED_DIC /   !<  H   a   AED2_CARBON_DATA_T%ID_FSED_CH4 +   i<  H   a   AED2_CARBON_DATA_T%ID_TEMP +   ±<  H   a   AED2_CARBON_DATA_T%ID_SALT +   ù<  H   a   AED2_CARBON_DATA_T%ID_WIND *   A=  H   a   AED2_CARBON_DATA_T%ID_VEL ,   =  H   a   AED2_CARBON_DATA_T%ID_DEPTH ,   Ñ=  H   a   AED2_CARBON_DATA_T%ID_CH4OX +   >  H   a   AED2_CARBON_DATA_T%ID_PCO2 .   a>  H   a   AED2_CARBON_DATA_T%ID_SED_DIC .   ©>  H   a   AED2_CARBON_DATA_T%ID_SED_CH4 2   ñ>  H   a   AED2_CARBON_DATA_T%ID_SED_CH4_EBB .   9?  H   a   AED2_CARBON_DATA_T%ID_ATM_CO2 .   ?  H   a   AED2_CARBON_DATA_T%ID_ATM_CH4 2   É?  H   a   AED2_CARBON_DATA_T%ID_ATM_CH4_EBB *   @  H   a   AED2_CARBON_DATA_T%ID_PAR +   Y@  H   a   AED2_CARBON_DATA_T%ID_EXTC )   ¡@  H   a   AED2_CARBON_DATA_T%ID_DZ *   é@  H   a   AED2_CARBON_DATA_T%ID_TAU ,   1A  H   a   AED2_CARBON_DATA_T%FSED_DIC ,   yA  H   a   AED2_CARBON_DATA_T%KSED_DIC 1   ÁA  H   a   AED2_CARBON_DATA_T%THETA_SED_DIC ,   	B  H   a   AED2_CARBON_DATA_T%FSED_CH4 ,   QB  H   a   AED2_CARBON_DATA_T%KSED_CH4 1   B  H   a   AED2_CARBON_DATA_T%THETA_SED_CH4 0   áB  H   a   AED2_CARBON_DATA_T%FSED_CH4_EBB 0   )C  H   a   AED2_CARBON_DATA_T%CH4_BUB_TAU0 *   qC  H   a   AED2_CARBON_DATA_T%RCH4OX *   ¹C  H   a   AED2_CARBON_DATA_T%KCH4OX +   D  H   a   AED2_CARBON_DATA_T%VTCH4OX +   ID  H   a   AED2_CARBON_DATA_T%ATM_CO2 +   D  H   a   AED2_CARBON_DATA_T%ATM_CH4 )   ÙD  H   a   AED2_CARBON_DATA_T%IONIC /   !E  H   a   AED2_CARBON_DATA_T%MAXMPBPRODN )   iE  H   a   AED2_CARBON_DATA_T%IKMPB +   ±E  H   a   AED2_CARBON_DATA_T%USE_OXY 5   ùE  H   a   AED2_CARBON_DATA_T%USE_SED_MODEL_DIC 5   AF  H   a   AED2_CARBON_DATA_T%USE_SED_MODEL_CH4 *   F  H   a   AED2_CARBON_DATA_T%SIMDIC *   ÑF  H   a   AED2_CARBON_DATA_T%SIMCH4 -   G  H   a   AED2_CARBON_DATA_T%SIMCH4EBB ,   aG  H   a   AED2_CARBON_DATA_T%ALK_MODE -   ©G  H   a   AED2_CARBON_DATA_T%CO2_MODEL 4   ñG  H   a   AED2_CARBON_DATA_T%CO2_PISTON_MODEL 4   9H  H   a   AED2_CARBON_DATA_T%CH4_PISTON_MODEL *   H  `   a   AED2_CARBON_DATA_T%DEFINE #   áH  ^      AED2_DEFINE_CARBON (   ?I  `   a   AED2_DEFINE_CARBON%DATA *   I  @   a   AED2_DEFINE_CARBON%NAMLST 5   ßI  k   a   AED2_CARBON_DATA_T%CALCULATE_SURFACE .   JJ  m      AED2_CALCULATE_SURFACE_CARBON 3   ·J  `   a   AED2_CALCULATE_SURFACE_CARBON%DATA 5   K     a   AED2_CALCULATE_SURFACE_CARBON%COLUMN 8   ¶K  @   a   AED2_CALCULATE_SURFACE_CARBON%LAYER_IDX -   öK  c   a   AED2_CARBON_DATA_T%CALCULATE &   YL  m      AED2_CALCULATE_CARBON +   ÆL  `   a   AED2_CALCULATE_CARBON%DATA -   &M     a   AED2_CALCULATE_CARBON%COLUMN 0   ÅM  @   a   AED2_CALCULATE_CARBON%LAYER_IDX 5   N  k   a   AED2_CARBON_DATA_T%CALCULATE_BENTHIC .   pN  m      AED2_CALCULATE_BENTHIC_CARBON 3   ÝN  `   a   AED2_CALCULATE_BENTHIC_CARBON%DATA 5   =O     a   AED2_CALCULATE_BENTHIC_CARBON%COLUMN 8   ÜO  @   a   AED2_CALCULATE_BENTHIC_CARBON%LAYER_IDX /   P  e   a   AED2_CARBON_DATA_T%EQUILIBRATE (   P  m      AED2_EQUILIBRATE_CARBON -   îP  `   a   AED2_EQUILIBRATE_CARBON%DATA /   NQ     a   AED2_EQUILIBRATE_CARBON%COLUMN 2   íQ  @   a   AED2_EQUILIBRATE_CARBON%LAYER_IDX 
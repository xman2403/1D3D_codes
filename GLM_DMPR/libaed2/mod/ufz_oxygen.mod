  ÉE  ´   k820309    ?          18.0        |#ó]                                                                                                          
       src/ufz_oxygen.F90 UFZ_OXYGEN              UFZ_OXYGEN_DATA_T                      @                              
                            @                              
       AED2_GAS_PISTON_VELOCITY AED2_OXYGEN_SAT                   @                               'È                    #AED2_MODEL_ID    #AED2_MODEL_NAME    #AED2_MODEL_PREFIX    #NEXT    #DEFINE    #INITIALIZE    #CALCULATE_SURFACE    #CALCULATE    #CALCULATE_BENTHIC    #CALCULATE_RIPARIAN !   #CALCULATE_DRY '   #EQUILIBRATE ,   #LIGHT_EXTINCTION 1   #RAIN_LOSS 7   #LIGHT_SHADING =   #BIO_DRAG C   #PARTICLE_BGC I   #MOBILITY P   #VALIDATE V   #DELETE [                 $                                                               $                                  @                                   $                                         D                          $                                  È       H             #AED2_MODEL_DATA_T                            ä              y#AED2_MODEL_DATA_T                                                   1         À    $                                              #AED2_DEFINE 	   #         @     @                            	                    #DATA 
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
                                |                    1 %         @                                }                           #NAME ~   #UNITS    #LONGNAME    #SURF              
                                ~                    1           
                                                    1           
                                                    1           
                                            %         @                                                           #NAME    #UNITS    #LONGNAME              
                                                    1           
                                                    1           
                                                    1 %         @                                                           #NAME              
                                                    1                                                  
                 
                                 0.                  @                              '8                   #AED2_MODEL_DATA_T    #ID_OXY    #ID_TEMP    #ID_SALT    #ID_DEPTH    #ID_LDEPTH    #ID_LAREA    #ID_LHT    #ID_WIND    #ID_FSED_OXY    #ID_OXY_SAT    #ID_ATM_OXY_EXCH    #ID_SED_OXY    #BOD_OXY    #SOD_OXY    #FSED_OXY    #KSED_OXY    #THETA_SED_OXY    #THETA_BOD_OXY    #THETA_SOD_OXY    #USE_LIMOD_OXYGEN_MODEL    #USE_SED_MODEL    #DEFINE     #CALCULATE_SURFACE ¤   #CALCULATE ©   #CALCULATE_BENTHIC ®                 $                                   È                      #AED2_MODEL_DATA_T                  $                                   È                           $                                   Ì                           $                                   Ð                           $                                   Ô                           $                                   Ø                           $                                   Ü                           $                                   à                           $                                   ä       	                    $                                   è       
                    $                                   ì                           $                                   ð                           $                                   ô                           $                                  ø          
                 $                                            
                 $                                           
                 $                                           
                 $                                           
                 $                                            
                 $                                  (         
                 $                                   0                          $                                   4            1         À    $                                              #AED2_DEFINE_OXYGEN ¡   #         @     @                             ¡                    #DATA ¢   #NAMLST £             
D                                ¢     8              #UFZ_OXYGEN_DATA_T              
                                  £           1         À    $                           ¤                  #AED2_CALCULATE_SURFACE_OXYGEN ¥   #         @     @                             ¥                    #DATA ¦   #COLUMN §   #LAYER_IDX ¨             
                                 ¦     8             #UFZ_OXYGEN_DATA_T              
D                                 §            °                       &                                           #AED2_COLUMN_T              
                                  ¨           1         À    $                           ©                  #AED2_CALCULATE_OXYGEN ª   #         @     @                             ª                    #DATA «   #COLUMN ¬   #LAYER_IDX ­             
                                 «     8             #UFZ_OXYGEN_DATA_T              
D                                 ¬            °                       &                                           #AED2_COLUMN_T              
                                  ­           1         À    $                           ®                  #AED2_CALCULATE_BENTHIC_OXYGEN ¯   #         @     @                             ¯                    #DATA °   #COLUMN ±   #LAYER_IDX ²             
                                 °     8             #UFZ_OXYGEN_DATA_T              
D                                 ±            °                       &                                           #AED2_COLUMN_T              
                                  ²                  &      fn#fn     Æ   "   b   uapp(UFZ_OXYGEN    è   @   J  AED2_CORE    (  i   J  AED2_UTIL ,     ®      AED2_MODEL_DATA_T+AED2_CORE :   ?  H   a   AED2_MODEL_DATA_T%AED2_MODEL_ID+AED2_CORE <     P   a   AED2_MODEL_DATA_T%AED2_MODEL_NAME+AED2_CORE >   ×  P   a   AED2_MODEL_DATA_T%AED2_MODEL_PREFIX+AED2_CORE 1   '  Þ   a   AED2_MODEL_DATA_T%NEXT+AED2_CORE 3     Y   a   AED2_MODEL_DATA_T%DEFINE+AED2_CORE &   ^  ^      AED2_DEFINE+AED2_CORE +   ¼  _   a   AED2_DEFINE%DATA+AED2_CORE -     @   a   AED2_DEFINE%NAMLST+AED2_CORE 7   [  ]   a   AED2_MODEL_DATA_T%INITIALIZE+AED2_CORE *   ¸  m      AED2_INITIALIZE+AED2_CORE /   %  _   a   AED2_INITIALIZE%DATA+AED2_CORE 1        a   AED2_INITIALIZE%COLUMN+AED2_CORE 4   #  @   a   AED2_INITIALIZE%LAYER_IDX+AED2_CORE >   c  d   a   AED2_MODEL_DATA_T%CALCULATE_SURFACE+AED2_CORE 1   Ç  m      AED2_CALCULATE_SURFACE+AED2_CORE 6   4	  _   a   AED2_CALCULATE_SURFACE%DATA+AED2_CORE 8   	     a   AED2_CALCULATE_SURFACE%COLUMN+AED2_CORE ;   2
  @   a   AED2_CALCULATE_SURFACE%LAYER_IDX+AED2_CORE 6   r
  \   a   AED2_MODEL_DATA_T%CALCULATE+AED2_CORE )   Î
  m      AED2_CALCULATE+AED2_CORE .   ;  _   a   AED2_CALCULATE%DATA+AED2_CORE 0        a   AED2_CALCULATE%COLUMN+AED2_CORE 3   9  @   a   AED2_CALCULATE%LAYER_IDX+AED2_CORE >   y  d   a   AED2_MODEL_DATA_T%CALCULATE_BENTHIC+AED2_CORE 1   Ý  m      AED2_CALCULATE_BENTHIC+AED2_CORE 6   J  _   a   AED2_CALCULATE_BENTHIC%DATA+AED2_CORE 8   ©     a   AED2_CALCULATE_BENTHIC%COLUMN+AED2_CORE ;   H  @   a   AED2_CALCULATE_BENTHIC%LAYER_IDX+AED2_CORE ?     e   a   AED2_MODEL_DATA_T%CALCULATE_RIPARIAN+AED2_CORE 2   í  y      AED2_CALCULATE_RIPARIAN+AED2_CORE 7   f  _   a   AED2_CALCULATE_RIPARIAN%DATA+AED2_CORE 9   Å     a   AED2_CALCULATE_RIPARIAN%COLUMN+AED2_CORE <   d  @   a   AED2_CALCULATE_RIPARIAN%LAYER_IDX+AED2_CORE 9   ¤  @   a   AED2_CALCULATE_RIPARIAN%PC_WET+AED2_CORE :   ä  `   a   AED2_MODEL_DATA_T%CALCULATE_DRY+AED2_CORE -   D  m      AED2_CALCULATE_DRY+AED2_CORE 2   ±  _   a   AED2_CALCULATE_DRY%DATA+AED2_CORE 4        a   AED2_CALCULATE_DRY%COLUMN+AED2_CORE 7   ¯  @   a   AED2_CALCULATE_DRY%LAYER_IDX+AED2_CORE 8   ï  ^   a   AED2_MODEL_DATA_T%EQUILIBRATE+AED2_CORE +   M  m      AED2_EQUILIBRATE+AED2_CORE 0   º  _   a   AED2_EQUILIBRATE%DATA+AED2_CORE 2        a   AED2_EQUILIBRATE%COLUMN+AED2_CORE 5   ¸  @   a   AED2_EQUILIBRATE%LAYER_IDX+AED2_CORE =   ø  c   a   AED2_MODEL_DATA_T%LIGHT_EXTINCTION+AED2_CORE 0   [  }      AED2_LIGHT_EXTINCTION+AED2_CORE 5   Ø  _   a   AED2_LIGHT_EXTINCTION%DATA+AED2_CORE 7   7     a   AED2_LIGHT_EXTINCTION%COLUMN+AED2_CORE :   Ö  @   a   AED2_LIGHT_EXTINCTION%LAYER_IDX+AED2_CORE ;     @   a   AED2_LIGHT_EXTINCTION%EXTINCTION+AED2_CORE 6   V  \   a   AED2_MODEL_DATA_T%RAIN_LOSS+AED2_CORE )   ²  x      AED2_RAIN_LOSS+AED2_CORE .   *  _   a   AED2_RAIN_LOSS%DATA+AED2_CORE 0        a   AED2_RAIN_LOSS%COLUMN+AED2_CORE 3   (  @   a   AED2_RAIN_LOSS%LAYER_IDX+AED2_CORE /   h  @   a   AED2_RAIN_LOSS%INFIL+AED2_CORE :   ¨  `   a   AED2_MODEL_DATA_T%LIGHT_SHADING+AED2_CORE -     }      AED2_LIGHT_SHADING+AED2_CORE 2     _   a   AED2_LIGHT_SHADING%DATA+AED2_CORE 4   ä     a   AED2_LIGHT_SHADING%COLUMN+AED2_CORE 7     @   a   AED2_LIGHT_SHADING%LAYER_IDX+AED2_CORE 8   Ã  @   a   AED2_LIGHT_SHADING%SHADE_FRAC+AED2_CORE 5     [   a   AED2_MODEL_DATA_T%BIO_DRAG+AED2_CORE (   ^  w      AED2_BIO_DRAG+AED2_CORE -   Õ  _   a   AED2_BIO_DRAG%DATA+AED2_CORE /   4     a   AED2_BIO_DRAG%COLUMN+AED2_CORE 2   Ó  @   a   AED2_BIO_DRAG%LAYER_IDX+AED2_CORE -     @   a   AED2_BIO_DRAG%DRAG+AED2_CORE 9   S  _   a   AED2_MODEL_DATA_T%PARTICLE_BGC+AED2_CORE ,   ²        AED2_PARTICLE_BGC+AED2_CORE 1   5  _   a   AED2_PARTICLE_BGC%DATA+AED2_CORE 3        a   AED2_PARTICLE_BGC%COLUMN+AED2_CORE 6   3   @   a   AED2_PARTICLE_BGC%LAYER_IDX+AED2_CORE 1   s   @   a   AED2_PARTICLE_BGC%PPID+AED2_CORE 3   ³      a   AED2_PARTICLE_BGC%PARTCL+AED2_CORE 5   ?!  [   a   AED2_MODEL_DATA_T%MOBILITY+AED2_CORE (   !  {      AED2_MOBILITY+AED2_CORE -   "  _   a   AED2_MOBILITY%DATA+AED2_CORE /   t"     a   AED2_MOBILITY%COLUMN+AED2_CORE 2   #  @   a   AED2_MOBILITY%LAYER_IDX+AED2_CORE 1   S#     a   AED2_MOBILITY%MOBILITY+AED2_CORE 5   ß#  [   a   AED2_MODEL_DATA_T%VALIDATE+AED2_CORE (   :$  u      AED2_VALIDATE+AED2_CORE -   ¯$  _   a   AED2_VALIDATE%DATA+AED2_CORE /   %     a   AED2_VALIDATE%COLUMN+AED2_CORE 2   ­%  @   a   AED2_VALIDATE%LAYER_IDX+AED2_CORE 3   í%  Y   a   AED2_MODEL_DATA_T%DELETE+AED2_CORE &   F&  R      AED2_DELETE+AED2_CORE +   &  _   a   AED2_DELETE%DATA+AED2_CORE (   ÷&  ¢       AED2_COLUMN_T+AED2_CORE -   '     a   AED2_COLUMN_T%CELL+AED2_CORE 3   -(  H   a   AED2_COLUMN_T%CELL_SHEET+AED2_CORE 1   u(  H   a   AED2_COLUMN_T%FLUX_ATM+AED2_CORE 1   ½(     a   AED2_COLUMN_T%FLUX_PEL+AED2_CORE 1   Q)  H   a   AED2_COLUMN_T%FLUX_BEN+AED2_CORE 1   )  H   a   AED2_COLUMN_T%FLUX_RIP+AED2_CORE 3   á)  ¸       AED2_GAS_PISTON_VELOCITY+AED2_UTIL 9   *  @   a   AED2_GAS_PISTON_VELOCITY%WSHGT+AED2_UTIL 8   Ù*  @   a   AED2_GAS_PISTON_VELOCITY%WIND+AED2_UTIL 7   +  @   a   AED2_GAS_PISTON_VELOCITY%TEM+AED2_UTIL 7   Y+  @   a   AED2_GAS_PISTON_VELOCITY%SAL+AED2_UTIL 7   +  @   a   AED2_GAS_PISTON_VELOCITY%VEL+AED2_UTIL 9   Ù+  @   a   AED2_GAS_PISTON_VELOCITY%DEPTH+AED2_UTIL 6   ,  @   a   AED2_GAS_PISTON_VELOCITY%LA+AED2_UTIL A   Y,  @   a   AED2_GAS_PISTON_VELOCITY%SCHMIDT_MODEL+AED2_UTIL @   ,  @   a   AED2_GAS_PISTON_VELOCITY%PISTON_MODEL+AED2_UTIL *   Ù,  d       AED2_OXYGEN_SAT+AED2_UTIL /   =-  @   a   AED2_OXYGEN_SAT%SALT+AED2_UTIL /   }-  @   a   AED2_OXYGEN_SAT%TEMP+AED2_UTIL    ½-  p       NAN_+AED2_CORE '   -.  v       SECS_PER_DAY+AED2_CORE /   £.  ¨       AED2_DEFINE_VARIABLE+AED2_CORE 4   K/  L   a   AED2_DEFINE_VARIABLE%NAME+AED2_CORE 5   /  L   a   AED2_DEFINE_VARIABLE%UNITS+AED2_CORE 8   ã/  L   a   AED2_DEFINE_VARIABLE%LONGNAME+AED2_CORE 7   /0  @   a   AED2_DEFINE_VARIABLE%INITIAL+AED2_CORE 7   o0  @   a   AED2_DEFINE_VARIABLE%MINIMUM+AED2_CORE 7   ¯0  @   a   AED2_DEFINE_VARIABLE%MAXIMUM+AED2_CORE 8   ï0  @   a   AED2_DEFINE_VARIABLE%MOBILITY+AED2_CORE 3   /1  Z       AED2_LOCATE_GLOBAL_SHEET+AED2_CORE 8   1  L   a   AED2_LOCATE_GLOBAL_SHEET%NAME+AED2_CORE :   Õ1  }       AED2_DEFINE_SHEET_DIAG_VARIABLE+AED2_CORE ?   R2  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%NAME+AED2_CORE @   2  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%UNITS+AED2_CORE C   ê2  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%LONGNAME+AED2_CORE ?   63  @   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%SURF+AED2_CORE 4   v3  s       AED2_DEFINE_DIAG_VARIABLE+AED2_CORE 9   é3  L   a   AED2_DEFINE_DIAG_VARIABLE%NAME+AED2_CORE :   54  L   a   AED2_DEFINE_DIAG_VARIABLE%UNITS+AED2_CORE =   4  L   a   AED2_DEFINE_DIAG_VARIABLE%LONGNAME+AED2_CORE -   Í4  Z       AED2_LOCATE_GLOBAL+AED2_CORE 2   '5  L   a   AED2_LOCATE_GLOBAL%NAME+AED2_CORE     s5  r       ZERO_+AED2_CORE "   å5  þ      UFZ_OXYGEN_DATA_T 4   ã7  g   a   UFZ_OXYGEN_DATA_T%AED2_MODEL_DATA_T )   J8  H   a   UFZ_OXYGEN_DATA_T%ID_OXY *   8  H   a   UFZ_OXYGEN_DATA_T%ID_TEMP *   Ú8  H   a   UFZ_OXYGEN_DATA_T%ID_SALT +   "9  H   a   UFZ_OXYGEN_DATA_T%ID_DEPTH ,   j9  H   a   UFZ_OXYGEN_DATA_T%ID_LDEPTH +   ²9  H   a   UFZ_OXYGEN_DATA_T%ID_LAREA )   ú9  H   a   UFZ_OXYGEN_DATA_T%ID_LHT *   B:  H   a   UFZ_OXYGEN_DATA_T%ID_WIND .   :  H   a   UFZ_OXYGEN_DATA_T%ID_FSED_OXY -   Ò:  H   a   UFZ_OXYGEN_DATA_T%ID_OXY_SAT 2   ;  H   a   UFZ_OXYGEN_DATA_T%ID_ATM_OXY_EXCH -   b;  H   a   UFZ_OXYGEN_DATA_T%ID_SED_OXY *   ª;  H   a   UFZ_OXYGEN_DATA_T%BOD_OXY *   ò;  H   a   UFZ_OXYGEN_DATA_T%SOD_OXY +   :<  H   a   UFZ_OXYGEN_DATA_T%FSED_OXY +   <  H   a   UFZ_OXYGEN_DATA_T%KSED_OXY 0   Ê<  H   a   UFZ_OXYGEN_DATA_T%THETA_SED_OXY 0   =  H   a   UFZ_OXYGEN_DATA_T%THETA_BOD_OXY 0   Z=  H   a   UFZ_OXYGEN_DATA_T%THETA_SOD_OXY 9   ¢=  H   a   UFZ_OXYGEN_DATA_T%USE_LIMOD_OXYGEN_MODEL 0   ê=  H   a   UFZ_OXYGEN_DATA_T%USE_SED_MODEL )   2>  `   a   UFZ_OXYGEN_DATA_T%DEFINE #   >  ^      AED2_DEFINE_OXYGEN (   ð>  _   a   AED2_DEFINE_OXYGEN%DATA *   O?  @   a   AED2_DEFINE_OXYGEN%NAMLST 4   ?  k   a   UFZ_OXYGEN_DATA_T%CALCULATE_SURFACE .   ú?  m      AED2_CALCULATE_SURFACE_OXYGEN 3   g@  _   a   AED2_CALCULATE_SURFACE_OXYGEN%DATA 5   Æ@     a   AED2_CALCULATE_SURFACE_OXYGEN%COLUMN 8   eA  @   a   AED2_CALCULATE_SURFACE_OXYGEN%LAYER_IDX ,   ¥A  c   a   UFZ_OXYGEN_DATA_T%CALCULATE &   B  m      AED2_CALCULATE_OXYGEN +   uB  _   a   AED2_CALCULATE_OXYGEN%DATA -   ÔB     a   AED2_CALCULATE_OXYGEN%COLUMN 0   sC  @   a   AED2_CALCULATE_OXYGEN%LAYER_IDX 4   ³C  k   a   UFZ_OXYGEN_DATA_T%CALCULATE_BENTHIC .   D  m      AED2_CALCULATE_BENTHIC_OXYGEN 3   D  _   a   AED2_CALCULATE_BENTHIC_OXYGEN%DATA 5   êD     a   AED2_CALCULATE_BENTHIC_OXYGEN%COLUMN 8   E  @   a   AED2_CALCULATE_BENTHIC_OXYGEN%LAYER_IDX 
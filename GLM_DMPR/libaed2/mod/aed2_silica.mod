  �:  �   k820309    ?          18.0        |#�]                                                                                                          
       src/aed2_silica.F90 AED2_SILICA              AED2_SILICA_DATA_T                      @                              
                         @                               '�                    #AED2_MODEL_ID    #AED2_MODEL_NAME    #AED2_MODEL_PREFIX    #NEXT    #DEFINE    #INITIALIZE    #CALCULATE_SURFACE    #CALCULATE    #CALCULATE_BENTHIC    #CALCULATE_RIPARIAN     #CALCULATE_DRY &   #EQUILIBRATE +   #LIGHT_EXTINCTION 0   #RAIN_LOSS 6   #LIGHT_SHADING <   #BIO_DRAG B   #PARTICLE_BGC H   #MOBILITY O   #VALIDATE U   #DELETE Z                � $                                                              � $                                  @                                  � $                                         D                          �$                                  �       H             #AED2_MODEL_DATA_T                            �              y#AED2_MODEL_DATA_T                                                   1         �   � $                      �                        #AED2_DEFINE    #         @     @                                                #DATA 	   #NAMLST 
             
                                	     �               #AED2_MODEL_DATA_T              
                                  
           1         �   � $                      �                        #AED2_INITIALIZE    #         @     @                                                #DATA    #COLUMN    #LAYER_IDX              
                                      �              #AED2_MODEL_DATA_T              
                                             �                       &                                           #AED2_COLUMN_T              
                                             1         �   � $                      �                        #AED2_CALCULATE_SURFACE    #         @     @                                                #DATA    #COLUMN    #LAYER_IDX              
                                      �              #AED2_MODEL_DATA_T              
                                             �                       &                                           #AED2_COLUMN_T              
                                             1         �   � $                      �                        #AED2_CALCULATE    #         @     @                                                #DATA    #COLUMN    #LAYER_IDX              
                                      �              #AED2_MODEL_DATA_T              
                                             �                       &                                           #AED2_COLUMN_T              
                                             1         �   � $                      �                   	     #AED2_CALCULATE_BENTHIC    #         @     @                                                #DATA    #COLUMN    #LAYER_IDX              
                                      �              #AED2_MODEL_DATA_T              
                                             �                       &                                           #AED2_COLUMN_T              
                                             1         �   � $                      �                    
     #AED2_CALCULATE_RIPARIAN !   #         @     @                            !                    #DATA "   #COLUMN #   #LAYER_IDX $   #PC_WET %             
                                 "     �              #AED2_MODEL_DATA_T              
                                 #            �       	                &                                           #AED2_COLUMN_T              
                                  $                     
                                 %     
      1         �   � $                      �      &                  #AED2_CALCULATE_DRY '   #         @     @                            '                    #DATA (   #COLUMN )   #LAYER_IDX *             
                                 (     �              #AED2_MODEL_DATA_T              
                                 )            �       
                &                                           #AED2_COLUMN_T              
                                  *           1         �   � $                      �      +                  #AED2_EQUILIBRATE ,   #         @     @                            ,                    #DATA -   #COLUMN .   #LAYER_IDX /             
                                 -     �              #AED2_MODEL_DATA_T              
                                 .            �                       &                                           #AED2_COLUMN_T              
                                  /           1         �   � $                      �      0              	    #AED2_LIGHT_EXTINCTION 1   #         @     @                            1                    #DATA 2   #COLUMN 3   #LAYER_IDX 4   #EXTINCTION 5             
                                 2     �              #AED2_MODEL_DATA_T              
                                 3            �                       &                                           #AED2_COLUMN_T              
                                  4                     
                                5     
       1         �   � $                      �      6              
    #AED2_RAIN_LOSS 7   #         @     @                            7                    #DATA 8   #COLUMN 9   #LAYER_IDX :   #INFIL ;             
                                 8     �              #AED2_MODEL_DATA_T              
                                 9            �                       &                                           #AED2_COLUMN_T              
                                  :                     
                                ;     
       1         �   � $                      �      <                  #AED2_LIGHT_SHADING =   #         @     @                            =                    #DATA >   #COLUMN ?   #LAYER_IDX @   #SHADE_FRAC A             
                                 >     �              #AED2_MODEL_DATA_T              
                                 ?            �                       &                                           #AED2_COLUMN_T              
                                  @                     
                                A     
       1         �   � $                      �      B                  #AED2_BIO_DRAG C   #         @     @                            C                    #DATA D   #COLUMN E   #LAYER_IDX F   #DRAG G             
                                 D     �              #AED2_MODEL_DATA_T              
                                 E            �                       &                                           #AED2_COLUMN_T              
                                  F                     
                                G     
       1         �   � $                      �      H                  #AED2_PARTICLE_BGC I   #         @     @                            I                    #DATA J   #COLUMN K   #LAYER_IDX L   #PPID M   #PARTCL N             
                                 J     �              #AED2_MODEL_DATA_T              
                                 K            �                       &                                           #AED2_COLUMN_T              
                                  L                     
                                 M                      
                                N                   
               &                                           1         �   � $                      �      O                  #AED2_MOBILITY P   #         @     @                            P                    #DATA Q   #COLUMN R   #LAYER_IDX S   #MOBILITY T             
                                 Q     �              #AED2_MODEL_DATA_T              
                                 R            �                       &                                           #AED2_COLUMN_T              
                                  S                     
                                T                   
               &                                           1         �   � $                     �      U                  #AED2_VALIDATE V   %         @   @                           V                           #DATA W   #COLUMN X   #LAYER_IDX Y             
                                 W     �              #AED2_MODEL_DATA_T              
                                 X            �                       &                                           #AED2_COLUMN_T              
                                  Y           1         �   � $                      �      Z                  #AED2_DELETE [   #         @     @                            [                    #DATA \             
                                \     �               #AED2_MODEL_DATA_T                      @                                '�                    #CELL ]   #CELL_SHEET ^   #FLUX_ATM _   #FLUX_PEL `   #FLUX_BEN a   #FLUX_RIP b               �$                             ]                              
            &                                                        �$                             ^     H          
                �$                             _     P          
               �$                             `            X                 
            &                                                        �$                             a     �          
                �$                             b     �          
                                               c     
                 
                                 0.                                            d     
                   
                        ��                                                    e     
                 
                      �@        86400.%         @                                f                           #NAME g   #UNITS h   #LONGNAME i   #INITIAL j   #MINIMUM k   #MAXIMUM l   #MOBILITY m             
                                g                    1           
                                h                    1           
                                i                    1           
                                j     
                
                                k     
                
                                l     
                
                                m     
      %         @                                n                           #NAME o             
                                o                    1 %         @                                p                           #NAME q             
                                q                    1 %         @                                r                           #NAME s   #UNITS t   #LONGNAME u   #SURF v             
                                s                    1           
                                t                    1           
                                u                    1           
                                 v           %         @                                w                           #NAME x             
                                x                    1                   @                         y     '                    #AED2_MODEL_DATA_T z   #ID_RSI {   #ID_OXY |   #ID_FSED_RSI }   #ID_TEMP ~   #ID_SED_RSI    #FSED_RSI �   #KSED_RSI �   #THETA_SED_RSI �   #USE_OXY �   #USE_SED_MODEL �   #DEFINE �   #CALCULATE �   #CALCULATE_BENTHIC �                � $                              z     �                      #AED2_MODEL_DATA_T                 � $                              {     �                          � $                              |     �                          � $                              }     �                          � $                              ~     �                          � $                                   �                          � $                             �     �          
                � $                             �     �          
                � $                             �     �       	   
                � $                              �     �       
                   � $                              �     �             1         �   � $                      �     �                  #AED2_DEFINE_SILICA �   #         @     @                             �                    #DATA �   #NAMLST �             
D                                �                    #AED2_SILICA_DATA_T y             
                                  �           1         �   � $                      �     �                  #AED2_CALCULATE_SILICA �   #         @     @                             �                    #DATA �   #COLUMN �   #LAYER_IDX �             
                                 �                   #AED2_SILICA_DATA_T y             
                                 �            �                       &                                           #AED2_COLUMN_T              
                                  �           1         �   � $                      �     �                  #AED2_CALCULATE_BENTHIC_SILICA �   #         @     @                             �                    #DATA �   #COLUMN �   #LAYER_IDX �             
                                 �                   #AED2_SILICA_DATA_T y             
D                                 �            �                       &                                           #AED2_COLUMN_T              
                                  �              �   (      fn#fn !   �   #   b   uapp(AED2_SILICA    �   @   J  AED2_CORE ,   +  �      AED2_MODEL_DATA_T+AED2_CORE :   �  H   a   AED2_MODEL_DATA_T%AED2_MODEL_ID+AED2_CORE <   !  P   a   AED2_MODEL_DATA_T%AED2_MODEL_NAME+AED2_CORE >   q  P   a   AED2_MODEL_DATA_T%AED2_MODEL_PREFIX+AED2_CORE 1   �  �   a   AED2_MODEL_DATA_T%NEXT+AED2_CORE 3   �  Y   a   AED2_MODEL_DATA_T%DEFINE+AED2_CORE &   �  ^      AED2_DEFINE+AED2_CORE +   V  _   a   AED2_DEFINE%DATA+AED2_CORE -   �  @   a   AED2_DEFINE%NAMLST+AED2_CORE 7   �  ]   a   AED2_MODEL_DATA_T%INITIALIZE+AED2_CORE *   R  m      AED2_INITIALIZE+AED2_CORE /   �  _   a   AED2_INITIALIZE%DATA+AED2_CORE 1     �   a   AED2_INITIALIZE%COLUMN+AED2_CORE 4   �  @   a   AED2_INITIALIZE%LAYER_IDX+AED2_CORE >   �  d   a   AED2_MODEL_DATA_T%CALCULATE_SURFACE+AED2_CORE 1   a  m      AED2_CALCULATE_SURFACE+AED2_CORE 6   �  _   a   AED2_CALCULATE_SURFACE%DATA+AED2_CORE 8   -	  �   a   AED2_CALCULATE_SURFACE%COLUMN+AED2_CORE ;   �	  @   a   AED2_CALCULATE_SURFACE%LAYER_IDX+AED2_CORE 6   
  \   a   AED2_MODEL_DATA_T%CALCULATE+AED2_CORE )   h
  m      AED2_CALCULATE+AED2_CORE .   �
  _   a   AED2_CALCULATE%DATA+AED2_CORE 0   4  �   a   AED2_CALCULATE%COLUMN+AED2_CORE 3   �  @   a   AED2_CALCULATE%LAYER_IDX+AED2_CORE >     d   a   AED2_MODEL_DATA_T%CALCULATE_BENTHIC+AED2_CORE 1   w  m      AED2_CALCULATE_BENTHIC+AED2_CORE 6   �  _   a   AED2_CALCULATE_BENTHIC%DATA+AED2_CORE 8   C  �   a   AED2_CALCULATE_BENTHIC%COLUMN+AED2_CORE ;   �  @   a   AED2_CALCULATE_BENTHIC%LAYER_IDX+AED2_CORE ?   "  e   a   AED2_MODEL_DATA_T%CALCULATE_RIPARIAN+AED2_CORE 2   �  y      AED2_CALCULATE_RIPARIAN+AED2_CORE 7      _   a   AED2_CALCULATE_RIPARIAN%DATA+AED2_CORE 9   _  �   a   AED2_CALCULATE_RIPARIAN%COLUMN+AED2_CORE <   �  @   a   AED2_CALCULATE_RIPARIAN%LAYER_IDX+AED2_CORE 9   >  @   a   AED2_CALCULATE_RIPARIAN%PC_WET+AED2_CORE :   ~  `   a   AED2_MODEL_DATA_T%CALCULATE_DRY+AED2_CORE -   �  m      AED2_CALCULATE_DRY+AED2_CORE 2   K  _   a   AED2_CALCULATE_DRY%DATA+AED2_CORE 4   �  �   a   AED2_CALCULATE_DRY%COLUMN+AED2_CORE 7   I  @   a   AED2_CALCULATE_DRY%LAYER_IDX+AED2_CORE 8   �  ^   a   AED2_MODEL_DATA_T%EQUILIBRATE+AED2_CORE +   �  m      AED2_EQUILIBRATE+AED2_CORE 0   T  _   a   AED2_EQUILIBRATE%DATA+AED2_CORE 2   �  �   a   AED2_EQUILIBRATE%COLUMN+AED2_CORE 5   R  @   a   AED2_EQUILIBRATE%LAYER_IDX+AED2_CORE =   �  c   a   AED2_MODEL_DATA_T%LIGHT_EXTINCTION+AED2_CORE 0   �  }      AED2_LIGHT_EXTINCTION+AED2_CORE 5   r  _   a   AED2_LIGHT_EXTINCTION%DATA+AED2_CORE 7   �  �   a   AED2_LIGHT_EXTINCTION%COLUMN+AED2_CORE :   p  @   a   AED2_LIGHT_EXTINCTION%LAYER_IDX+AED2_CORE ;   �  @   a   AED2_LIGHT_EXTINCTION%EXTINCTION+AED2_CORE 6   �  \   a   AED2_MODEL_DATA_T%RAIN_LOSS+AED2_CORE )   L  x      AED2_RAIN_LOSS+AED2_CORE .   �  _   a   AED2_RAIN_LOSS%DATA+AED2_CORE 0   #  �   a   AED2_RAIN_LOSS%COLUMN+AED2_CORE 3   �  @   a   AED2_RAIN_LOSS%LAYER_IDX+AED2_CORE /     @   a   AED2_RAIN_LOSS%INFIL+AED2_CORE :   B  `   a   AED2_MODEL_DATA_T%LIGHT_SHADING+AED2_CORE -   �  }      AED2_LIGHT_SHADING+AED2_CORE 2     _   a   AED2_LIGHT_SHADING%DATA+AED2_CORE 4   ~  �   a   AED2_LIGHT_SHADING%COLUMN+AED2_CORE 7     @   a   AED2_LIGHT_SHADING%LAYER_IDX+AED2_CORE 8   ]  @   a   AED2_LIGHT_SHADING%SHADE_FRAC+AED2_CORE 5   �  [   a   AED2_MODEL_DATA_T%BIO_DRAG+AED2_CORE (   �  w      AED2_BIO_DRAG+AED2_CORE -   o  _   a   AED2_BIO_DRAG%DATA+AED2_CORE /   �  �   a   AED2_BIO_DRAG%COLUMN+AED2_CORE 2   m  @   a   AED2_BIO_DRAG%LAYER_IDX+AED2_CORE -   �  @   a   AED2_BIO_DRAG%DRAG+AED2_CORE 9   �  _   a   AED2_MODEL_DATA_T%PARTICLE_BGC+AED2_CORE ,   L  �      AED2_PARTICLE_BGC+AED2_CORE 1   �  _   a   AED2_PARTICLE_BGC%DATA+AED2_CORE 3   .  �   a   AED2_PARTICLE_BGC%COLUMN+AED2_CORE 6   �  @   a   AED2_PARTICLE_BGC%LAYER_IDX+AED2_CORE 1      @   a   AED2_PARTICLE_BGC%PPID+AED2_CORE 3   M   �   a   AED2_PARTICLE_BGC%PARTCL+AED2_CORE 5   �   [   a   AED2_MODEL_DATA_T%MOBILITY+AED2_CORE (   4!  {      AED2_MOBILITY+AED2_CORE -   �!  _   a   AED2_MOBILITY%DATA+AED2_CORE /   "  �   a   AED2_MOBILITY%COLUMN+AED2_CORE 2   �"  @   a   AED2_MOBILITY%LAYER_IDX+AED2_CORE 1   �"  �   a   AED2_MOBILITY%MOBILITY+AED2_CORE 5   y#  [   a   AED2_MODEL_DATA_T%VALIDATE+AED2_CORE (   �#  u      AED2_VALIDATE+AED2_CORE -   I$  _   a   AED2_VALIDATE%DATA+AED2_CORE /   �$  �   a   AED2_VALIDATE%COLUMN+AED2_CORE 2   G%  @   a   AED2_VALIDATE%LAYER_IDX+AED2_CORE 3   �%  Y   a   AED2_MODEL_DATA_T%DELETE+AED2_CORE &   �%  R      AED2_DELETE+AED2_CORE +   2&  _   a   AED2_DELETE%DATA+AED2_CORE (   �&  �       AED2_COLUMN_T+AED2_CORE -   3'  �   a   AED2_COLUMN_T%CELL+AED2_CORE 3   �'  H   a   AED2_COLUMN_T%CELL_SHEET+AED2_CORE 1   (  H   a   AED2_COLUMN_T%FLUX_ATM+AED2_CORE 1   W(  �   a   AED2_COLUMN_T%FLUX_PEL+AED2_CORE 1   �(  H   a   AED2_COLUMN_T%FLUX_BEN+AED2_CORE 1   3)  H   a   AED2_COLUMN_T%FLUX_RIP+AED2_CORE     {)  r       ZERO_+AED2_CORE    �)  p       NAN_+AED2_CORE '   ]*  v       SECS_PER_DAY+AED2_CORE /   �*  �       AED2_DEFINE_VARIABLE+AED2_CORE 4   {+  L   a   AED2_DEFINE_VARIABLE%NAME+AED2_CORE 5   �+  L   a   AED2_DEFINE_VARIABLE%UNITS+AED2_CORE 8   ,  L   a   AED2_DEFINE_VARIABLE%LONGNAME+AED2_CORE 7   _,  @   a   AED2_DEFINE_VARIABLE%INITIAL+AED2_CORE 7   �,  @   a   AED2_DEFINE_VARIABLE%MINIMUM+AED2_CORE 7   �,  @   a   AED2_DEFINE_VARIABLE%MAXIMUM+AED2_CORE 8   -  @   a   AED2_DEFINE_VARIABLE%MOBILITY+AED2_CORE /   _-  Z       AED2_LOCATE_VARIABLE+AED2_CORE 4   �-  L   a   AED2_LOCATE_VARIABLE%NAME+AED2_CORE 3   .  Z       AED2_LOCATE_GLOBAL_SHEET+AED2_CORE 8   _.  L   a   AED2_LOCATE_GLOBAL_SHEET%NAME+AED2_CORE :   �.  }       AED2_DEFINE_SHEET_DIAG_VARIABLE+AED2_CORE ?   (/  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%NAME+AED2_CORE @   t/  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%UNITS+AED2_CORE C   �/  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%LONGNAME+AED2_CORE ?   0  @   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%SURF+AED2_CORE -   L0  Z       AED2_LOCATE_GLOBAL+AED2_CORE 2   �0  L   a   AED2_LOCATE_GLOBAL%NAME+AED2_CORE #   �0  .      AED2_SILICA_DATA_T 5    2  g   a   AED2_SILICA_DATA_T%AED2_MODEL_DATA_T *   �2  H   a   AED2_SILICA_DATA_T%ID_RSI *   �2  H   a   AED2_SILICA_DATA_T%ID_OXY /   3  H   a   AED2_SILICA_DATA_T%ID_FSED_RSI +   _3  H   a   AED2_SILICA_DATA_T%ID_TEMP .   �3  H   a   AED2_SILICA_DATA_T%ID_SED_RSI ,   �3  H   a   AED2_SILICA_DATA_T%FSED_RSI ,   74  H   a   AED2_SILICA_DATA_T%KSED_RSI 1   4  H   a   AED2_SILICA_DATA_T%THETA_SED_RSI +   �4  H   a   AED2_SILICA_DATA_T%USE_OXY 1   5  H   a   AED2_SILICA_DATA_T%USE_SED_MODEL *   W5  `   a   AED2_SILICA_DATA_T%DEFINE #   �5  ^      AED2_DEFINE_SILICA (   6  `   a   AED2_DEFINE_SILICA%DATA *   u6  @   a   AED2_DEFINE_SILICA%NAMLST -   �6  c   a   AED2_SILICA_DATA_T%CALCULATE &   7  m      AED2_CALCULATE_SILICA +   �7  `   a   AED2_CALCULATE_SILICA%DATA -   �7  �   a   AED2_CALCULATE_SILICA%COLUMN 0   �8  @   a   AED2_CALCULATE_SILICA%LAYER_IDX 5   �8  k   a   AED2_SILICA_DATA_T%CALCULATE_BENTHIC .   /9  m      AED2_CALCULATE_BENTHIC_SILICA 3   �9  `   a   AED2_CALCULATE_BENTHIC_SILICA%DATA 5   �9  �   a   AED2_CALCULATE_BENTHIC_SILICA%COLUMN 8   �:  @   a   AED2_CALCULATE_BENTHIC_SILICA%LAYER_IDX 
  uA  �   k820309    ?          18.0        {#�]                                                                                                          
       src/aed2_core.F90 AED2_CORE              AED2_MODEL_DATA_T AED2_VARIABLE_T AED2_COLUMN_T AED2_INIT_CORE AED2_GET_VAR AED2_CORE_STATUS AED2_SET_PREFIX AED2_DEFINE_VARIABLE AED2_DEFINE_SHEET_VARIABLE AED2_LOCATE_SHEET_VARIABLE AED2_LOCATE_VARIABLE AED2_DEFINE_DIAG_VARIABLE AED2_DEFINE_SHEET_DIAG_VARIABLE AED2_LOCATE_GLOBAL AED2_LOCATE_GLOBAL_SHEET HOST_HAS_CELL_VEL ZERO_ ONE_ NAN_ MISVAL_ SECS_PER_DAY CUR_MODEL_NAME                   @                               '�                    #AED2_MODEL_ID    #AED2_MODEL_NAME    #AED2_MODEL_PREFIX    #NEXT    #DEFINE    #INITIALIZE 
   #CALCULATE_SURFACE    #CALCULATE    #CALCULATE_BENTHIC    #CALCULATE_RIPARIAN    #CALCULATE_DRY %   #EQUILIBRATE *   #LIGHT_EXTINCTION /   #RAIN_LOSS 5   #LIGHT_SHADING ;   #BIO_DRAG A   #PARTICLE_BGC G   #MOBILITY N   #VALIDATE T   #DELETE Y                � $                                                              � $                                  @                                  � $                                         D                          �$                                  �       H             #AED2_MODEL_DATA_T                            �              y#AED2_MODEL_DATA_T                                                   1         �   � $                      �                        #AED2_DEFINE    #         @     @                                                 #DATA    #NAMLST 	             
                                     �               #AED2_MODEL_DATA_T              
                                  	           1         �   � $                      �      
                  #AED2_INITIALIZE    #         @     @                                                 #DATA    #COLUMN    #LAYER_IDX              
                                      �              #AED2_MODEL_DATA_T              
                                             �                       &                                           #AED2_COLUMN_T              
                                             1         �   � $                      �                        #AED2_CALCULATE_SURFACE    #         @     @                                                 #DATA    #COLUMN    #LAYER_IDX              
                                      �              #AED2_MODEL_DATA_T              
                                             �                       &                                           #AED2_COLUMN_T              
                                             1         �   � $                      �                        #AED2_CALCULATE    #         @     @                                                 #DATA    #COLUMN    #LAYER_IDX              
                                      �              #AED2_MODEL_DATA_T              
                                             �                       &                                           #AED2_COLUMN_T              
                                             1         �   � $                      �                   	     #AED2_CALCULATE_BENTHIC    #         @     @                                                 #DATA    #COLUMN    #LAYER_IDX              
                                      �              #AED2_MODEL_DATA_T              
                                             �                       &                                           #AED2_COLUMN_T              
                                             1         �   � $                      �                   
     #AED2_CALCULATE_RIPARIAN     #         @     @                                                  #DATA !   #COLUMN "   #LAYER_IDX #   #PC_WET $             
                                 !     �              #AED2_MODEL_DATA_T              
                                 "            �       	                &                                           #AED2_COLUMN_T              
                                  #                     
                                 $     
      1         �   � $                      �      %                  #AED2_CALCULATE_DRY &   #         @     @                             &                    #DATA '   #COLUMN (   #LAYER_IDX )             
                                 '     �              #AED2_MODEL_DATA_T              
                                 (            �       
                &                                           #AED2_COLUMN_T              
                                  )           1         �   � $                      �      *                  #AED2_EQUILIBRATE +   #         @     @                             +                    #DATA ,   #COLUMN -   #LAYER_IDX .             
                                 ,     �              #AED2_MODEL_DATA_T              
                                 -            �                       &                                           #AED2_COLUMN_T              
                                  .           1         �   � $                      �      /              	    #AED2_LIGHT_EXTINCTION 0   #         @     @                             0                    #DATA 1   #COLUMN 2   #LAYER_IDX 3   #EXTINCTION 4             
                                 1     �              #AED2_MODEL_DATA_T              
                                 2            �                       &                                           #AED2_COLUMN_T              
                                  3                     
                                4     
       1         �   � $                      �      5              
    #AED2_RAIN_LOSS 6   #         @     @                             6                    #DATA 7   #COLUMN 8   #LAYER_IDX 9   #INFIL :             
                                 7     �              #AED2_MODEL_DATA_T              
                                 8            �                       &                                           #AED2_COLUMN_T              
                                  9                     
                                :     
       1         �   � $                      �      ;                  #AED2_LIGHT_SHADING <   #         @     @                             <                    #DATA =   #COLUMN >   #LAYER_IDX ?   #SHADE_FRAC @             
                                 =     �              #AED2_MODEL_DATA_T              
                                 >            �                       &                                           #AED2_COLUMN_T              
                                  ?                     
                                @     
       1         �   � $                      �      A                  #AED2_BIO_DRAG B   #         @     @                             B                    #DATA C   #COLUMN D   #LAYER_IDX E   #DRAG F             
                                 C     �              #AED2_MODEL_DATA_T              
                                 D            �                       &                                           #AED2_COLUMN_T              
                                  E                     
                                F     
       1         �   � $                      �      G                  #AED2_PARTICLE_BGC H   #         @     @                             H                    #DATA I   #COLUMN J   #LAYER_IDX K   #PPID L   #PARTCL M             
                                 I     �              #AED2_MODEL_DATA_T              
                                 J            �                       &                                           #AED2_COLUMN_T              
                                  K                     
                                 L                      
                                M                   
               &                                           1         �   � $                      �      N                  #AED2_MOBILITY O   #         @     @                             O                    #DATA P   #COLUMN Q   #LAYER_IDX R   #MOBILITY S             
                                 P     �              #AED2_MODEL_DATA_T              
                                 Q            �                       &                                           #AED2_COLUMN_T              
                                  R                     
                                S                   
               &                                           1         �   � $                     �      T                  #AED2_VALIDATE U   %         @   @                            U                           #DATA V   #COLUMN W   #LAYER_IDX X             
                                 V     �              #AED2_MODEL_DATA_T              
                                 W            �                       &                                           #AED2_COLUMN_T              
                                  X           1         �   � $                      �      Y                  #AED2_DELETE Z   #         @     @                             Z                    #DATA [             
                                [     �               #AED2_MODEL_DATA_T                      @                           \     'X                   #NAME ]   #MODEL_NAME ^   #LONGNAME _   #UNITS `   #INITIAL a   #MINIMUM b   #MAXIMUM c   #MOBILITY d   #LIGHT_EXTINCTION e   #SHEET f   #DIAG g   #EXTERN h   #FOUND i   #TOP j   #BOT k                � $                             ]     @                                   � $                             ^     @       @                           � $                             _     �       �                           � $                             `                                       � $                             a              
                � $                             b               
                � $                             c     (         
                � $                             d     0         
                � $                             e     8      	   
                � $                              f     @      
                   � $                              g     D                         � $                              h     H                         � $                              i     L                         � $                              j     P                         � $                              k     T                              @                                '�                    #CELL l   #CELL_SHEET m   #FLUX_ATM n   #FLUX_PEL o   #FLUX_BEN p   #FLUX_RIP q               �$                             l                              
            &                                                        �$                             m     H          
                �$                             n     P          
               �$                             o            X                 
            &                                                        �$                             p     �          
                �$                             q     �          
   %         @                                r                           #DNAME s   #HAVE_CELL_VEL t                                            s                     1            @                               t            %         @                                u                           #WHICH v   #TVAR w             
                                  v                     D                                w     X              #AED2_VARIABLE_T \   %         @                                x                           #N_V y   #N_SV z   #N_D {   #N_SD |             D                                 y                      D                                 z                      D                                 {                      D                                 |            #         @                                   }                    #PREFIX ~                                            ~                            %         @                                                            #NAME �   #UNITS �   #LONGNAME �   #INITIAL �   #MINIMUM �   #MAXIMUM �   #MOBILITY �             
  @                             �                    1           
  @                             �                    1           
  @                             �                    1           
 @                              �     
                
 @                              �     
                
 @                              �     
                
 @                              �     
      %         @                                 �                           #NAME �   #UNITS �   #LONGNAME �   #INITIAL �   #MINIMUM �   #MAXIMUM �   #SURF �             
  @                             �                    1           
  @                             �                    1           
  @                             �                    1           
 @                              �     
                
 @                              �     
                
 @                              �     
                
 @                               �           %         @                                 �                           #NAME �             
  @                             �                    1 %         @                                 �                           #NAME �             
  @                             �                    1 %         @                                 �                           #NAME �   #UNITS �   #LONGNAME �             
  @                             �                    1           
  @                             �                    1           
  @                             �                    1 %         @                                 �                           #NAME �   #UNITS �   #LONGNAME �   #SURF �             
  @                             �                    1           
  @                             �                    1           
  @                             �                    1           
 @                               �           %         @                                 �                           #NAME �             
  @                             �                    1 %         @                                 �                           #NAME �             
  @                             �                    1           @                                 �                                                        �     
                 
                                 0.                                            �     
                 
                       �?        1.                                            �     
                   
                        ��                                                    �     
                   
                      ����                                                    �     
                 
                      �@        86400.                                           �     @          �   $      fn#fn    �   �  b   uapp(AED2_CORE "   M  �      AED2_MODEL_DATA_T 0   �  H   a   AED2_MODEL_DATA_T%AED2_MODEL_ID 2   C  P   a   AED2_MODEL_DATA_T%AED2_MODEL_NAME 4   �  P   a   AED2_MODEL_DATA_T%AED2_MODEL_PREFIX '   �  �   a   AED2_MODEL_DATA_T%NEXT )   �  Y   a   AED2_MODEL_DATA_T%DEFINE      ^      AED2_DEFINE !   x  _   a   AED2_DEFINE%DATA #   �  @   a   AED2_DEFINE%NAMLST -     ]   a   AED2_MODEL_DATA_T%INITIALIZE     t  m      AED2_INITIALIZE %   �  _   a   AED2_INITIALIZE%DATA '   @  �   a   AED2_INITIALIZE%COLUMN *   �  @   a   AED2_INITIALIZE%LAYER_IDX 4   	  d   a   AED2_MODEL_DATA_T%CALCULATE_SURFACE '   �	  m      AED2_CALCULATE_SURFACE ,   �	  _   a   AED2_CALCULATE_SURFACE%DATA .   O
  �   a   AED2_CALCULATE_SURFACE%COLUMN 1   �
  @   a   AED2_CALCULATE_SURFACE%LAYER_IDX ,   .  \   a   AED2_MODEL_DATA_T%CALCULATE    �  m      AED2_CALCULATE $   �  _   a   AED2_CALCULATE%DATA &   V  �   a   AED2_CALCULATE%COLUMN )   �  @   a   AED2_CALCULATE%LAYER_IDX 4   5  d   a   AED2_MODEL_DATA_T%CALCULATE_BENTHIC '   �  m      AED2_CALCULATE_BENTHIC ,     _   a   AED2_CALCULATE_BENTHIC%DATA .   e  �   a   AED2_CALCULATE_BENTHIC%COLUMN 1     @   a   AED2_CALCULATE_BENTHIC%LAYER_IDX 5   D  e   a   AED2_MODEL_DATA_T%CALCULATE_RIPARIAN (   �  y      AED2_CALCULATE_RIPARIAN -   "  _   a   AED2_CALCULATE_RIPARIAN%DATA /   �  �   a   AED2_CALCULATE_RIPARIAN%COLUMN 2      @   a   AED2_CALCULATE_RIPARIAN%LAYER_IDX /   `  @   a   AED2_CALCULATE_RIPARIAN%PC_WET 0   �  `   a   AED2_MODEL_DATA_T%CALCULATE_DRY #      m      AED2_CALCULATE_DRY (   m  _   a   AED2_CALCULATE_DRY%DATA *   �  �   a   AED2_CALCULATE_DRY%COLUMN -   k  @   a   AED2_CALCULATE_DRY%LAYER_IDX .   �  ^   a   AED2_MODEL_DATA_T%EQUILIBRATE !   	  m      AED2_EQUILIBRATE &   v  _   a   AED2_EQUILIBRATE%DATA (   �  �   a   AED2_EQUILIBRATE%COLUMN +   t  @   a   AED2_EQUILIBRATE%LAYER_IDX 3   �  c   a   AED2_MODEL_DATA_T%LIGHT_EXTINCTION &     }      AED2_LIGHT_EXTINCTION +   �  _   a   AED2_LIGHT_EXTINCTION%DATA -   �  �   a   AED2_LIGHT_EXTINCTION%COLUMN 0   �  @   a   AED2_LIGHT_EXTINCTION%LAYER_IDX 1   �  @   a   AED2_LIGHT_EXTINCTION%EXTINCTION ,     \   a   AED2_MODEL_DATA_T%RAIN_LOSS    n  x      AED2_RAIN_LOSS $   �  _   a   AED2_RAIN_LOSS%DATA &   E  �   a   AED2_RAIN_LOSS%COLUMN )   �  @   a   AED2_RAIN_LOSS%LAYER_IDX %   $  @   a   AED2_RAIN_LOSS%INFIL 0   d  `   a   AED2_MODEL_DATA_T%LIGHT_SHADING #   �  }      AED2_LIGHT_SHADING (   A  _   a   AED2_LIGHT_SHADING%DATA *   �  �   a   AED2_LIGHT_SHADING%COLUMN -   ?  @   a   AED2_LIGHT_SHADING%LAYER_IDX .     @   a   AED2_LIGHT_SHADING%SHADE_FRAC +   �  [   a   AED2_MODEL_DATA_T%BIO_DRAG      w      AED2_BIO_DRAG #   �  _   a   AED2_BIO_DRAG%DATA %   �  �   a   AED2_BIO_DRAG%COLUMN (   �  @   a   AED2_BIO_DRAG%LAYER_IDX #   �  @   a   AED2_BIO_DRAG%DRAG /     _   a   AED2_MODEL_DATA_T%PARTICLE_BGC "   n  �      AED2_PARTICLE_BGC '   �  _   a   AED2_PARTICLE_BGC%DATA )   P   �   a   AED2_PARTICLE_BGC%COLUMN ,   �   @   a   AED2_PARTICLE_BGC%LAYER_IDX '   /!  @   a   AED2_PARTICLE_BGC%PPID )   o!  �   a   AED2_PARTICLE_BGC%PARTCL +   �!  [   a   AED2_MODEL_DATA_T%MOBILITY    V"  {      AED2_MOBILITY #   �"  _   a   AED2_MOBILITY%DATA %   0#  �   a   AED2_MOBILITY%COLUMN (   �#  @   a   AED2_MOBILITY%LAYER_IDX '   $  �   a   AED2_MOBILITY%MOBILITY +   �$  [   a   AED2_MODEL_DATA_T%VALIDATE    �$  u      AED2_VALIDATE #   k%  _   a   AED2_VALIDATE%DATA %   �%  �   a   AED2_VALIDATE%COLUMN (   i&  @   a   AED2_VALIDATE%LAYER_IDX )   �&  Y   a   AED2_MODEL_DATA_T%DELETE    '  R      AED2_DELETE !   T'  _   a   AED2_DELETE%DATA     �'        AED2_VARIABLE_T %   �(  P   a   AED2_VARIABLE_T%NAME +   )  P   a   AED2_VARIABLE_T%MODEL_NAME )   _)  P   a   AED2_VARIABLE_T%LONGNAME &   �)  P   a   AED2_VARIABLE_T%UNITS (   �)  H   a   AED2_VARIABLE_T%INITIAL (   G*  H   a   AED2_VARIABLE_T%MINIMUM (   �*  H   a   AED2_VARIABLE_T%MAXIMUM )   �*  H   a   AED2_VARIABLE_T%MOBILITY 1   +  H   a   AED2_VARIABLE_T%LIGHT_EXTINCTION &   g+  H   a   AED2_VARIABLE_T%SHEET %   �+  H   a   AED2_VARIABLE_T%DIAG '   �+  H   a   AED2_VARIABLE_T%EXTERN &   ?,  H   a   AED2_VARIABLE_T%FOUND $   �,  H   a   AED2_VARIABLE_T%TOP $   �,  H   a   AED2_VARIABLE_T%BOT    -  �       AED2_COLUMN_T #   �-  �   a   AED2_COLUMN_T%CELL )   M.  H   a   AED2_COLUMN_T%CELL_SHEET '   �.  H   a   AED2_COLUMN_T%FLUX_ATM '   �.  �   a   AED2_COLUMN_T%FLUX_PEL '   q/  H   a   AED2_COLUMN_T%FLUX_BEN '   �/  H   a   AED2_COLUMN_T%FLUX_RIP    0  n       AED2_INIT_CORE %   o0  L   a   AED2_INIT_CORE%DNAME -   �0  @   a   AED2_INIT_CORE%HAVE_CELL_VEL    �0  e       AED2_GET_VAR #   `1  @   a   AED2_GET_VAR%WHICH "   �1  ]   a   AED2_GET_VAR%TVAR !   �1  v       AED2_CORE_STATUS %   s2  @   a   AED2_CORE_STATUS%N_V &   �2  @   a   AED2_CORE_STATUS%N_SV %   �2  @   a   AED2_CORE_STATUS%N_D &   33  @   a   AED2_CORE_STATUS%N_SD     s3  T       AED2_SET_PREFIX '   �3  P   a   AED2_SET_PREFIX%PREFIX %   4  �       AED2_DEFINE_VARIABLE *   �4  L   a   AED2_DEFINE_VARIABLE%NAME +   5  L   a   AED2_DEFINE_VARIABLE%UNITS .   W5  L   a   AED2_DEFINE_VARIABLE%LONGNAME -   �5  @   a   AED2_DEFINE_VARIABLE%INITIAL -   �5  @   a   AED2_DEFINE_VARIABLE%MINIMUM -   #6  @   a   AED2_DEFINE_VARIABLE%MAXIMUM .   c6  @   a   AED2_DEFINE_VARIABLE%MOBILITY +   �6  �       AED2_DEFINE_SHEET_VARIABLE 0   G7  L   a   AED2_DEFINE_SHEET_VARIABLE%NAME 1   �7  L   a   AED2_DEFINE_SHEET_VARIABLE%UNITS 4   �7  L   a   AED2_DEFINE_SHEET_VARIABLE%LONGNAME 3   +8  @   a   AED2_DEFINE_SHEET_VARIABLE%INITIAL 3   k8  @   a   AED2_DEFINE_SHEET_VARIABLE%MINIMUM 3   �8  @   a   AED2_DEFINE_SHEET_VARIABLE%MAXIMUM 0   �8  @   a   AED2_DEFINE_SHEET_VARIABLE%SURF +   +9  Z       AED2_LOCATE_SHEET_VARIABLE 0   �9  L   a   AED2_LOCATE_SHEET_VARIABLE%NAME %   �9  Z       AED2_LOCATE_VARIABLE *   +:  L   a   AED2_LOCATE_VARIABLE%NAME *   w:  s       AED2_DEFINE_DIAG_VARIABLE /   �:  L   a   AED2_DEFINE_DIAG_VARIABLE%NAME 0   6;  L   a   AED2_DEFINE_DIAG_VARIABLE%UNITS 3   �;  L   a   AED2_DEFINE_DIAG_VARIABLE%LONGNAME 0   �;  }       AED2_DEFINE_SHEET_DIAG_VARIABLE 5   K<  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%NAME 6   �<  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%UNITS 9   �<  L   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%LONGNAME 5   /=  @   a   AED2_DEFINE_SHEET_DIAG_VARIABLE%SURF #   o=  Z       AED2_LOCATE_GLOBAL (   �=  L   a   AED2_LOCATE_GLOBAL%NAME )   >  Z       AED2_LOCATE_GLOBAL_SHEET .   o>  L   a   AED2_LOCATE_GLOBAL_SHEET%NAME "   �>  @       HOST_HAS_CELL_VEL    �>  r       ZERO_    m?  r       ONE_    �?  p       NAN_    O@  p       MISVAL_    �@  v       SECS_PER_DAY    5A  @       CUR_MODEL_NAME 
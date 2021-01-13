                #Interpolate velocity
                #INCREMENT GLOBAL VELOCITY FROM INCREMENTS!!!
                #CHANGE F TO ASSEMBLE IN THE SAME PLACE FOR BOTH FORMS
                # THIS IS THE CRITERIA IN WHICH VARS ARE INBLOCK PER NODE
                # juf=0
                # uvf=0
                # for n in range (4):
                    # d=elnodes.astype(int)[e][n]
                    # for i in range (var_dim[0]):    #Velocity is var 0
                        # print("UV loc glob ",i,ndof*d+i)
                        # UV[i,0]=Uglob[ndof*d+i]
                    # uvf+=var_dim[0]
                    # for j in range (var_dim[1]):
                        # #print("J",j)
                        # if (form==1):
                            # Usig[j+juf,0]=Uglob[ndof*d+var_dim[0]+j]
                            # UF  [j+juf,0]=Uglob[ndof*d+6+j]
                        # else: #Fig 4.1, Z is not translated to Fvpt
                            # UF  [j+juf,0]=Uglob[ndof*d+var_dim[0]+j]
                            # #print("UF(j,coord)",j,ndof*d+6+j)
                    # juf+=var_dim[1]
                    
                    # if plastic:
                        # for j in range (var_dim[2]):
                            # UFvp[j,0]=Uglob[ndof*d+var_dim[0]+var_dim[1]+j]
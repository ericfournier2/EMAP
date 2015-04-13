bash Scripts/FullAnalysis.sh "Rachele0.1/Epigenetic" "Epigenetic.target" "Raw/Rachele" \
                     "Rachele0.1/Transcriptomic" "Transcriptomic.target" "Raw/Rachele" \
                     "Control" "Rachele0.1/Combined" 1.4 0.1 FALSE 1.5 0.05 FALSE v2

bash Scripts/FullAnalysis.sh "Rachele/V1/Epigenetic" "Design.target" "Raw/Rachele/V1" \
                     "" "" "" \
                     "Control" "" 1.5 0.05 FALSE 1.5 0.05 FALSE v1
                     
bash Scripts/FullAnalysis.sh "Rachele/V2/Epigenetic" "Design.target" "Raw/Rachele/V2" \
                     "" "" "" \
                     "Control" "" 1.5 0.05 FALSE 1.5 0.05 FALSE v2                     
                     
bash Scripts/FullAnalysis.sh "Rachele/V2/Epigenetic" "Design-All.target" "Raw/Rachele/V2" \
                     "" "" "" \
                     "Control" "" 1.5 0.05 FALSE 1.5 0.05 FALSE v2                                          
                     
bash Scripts/FullAnalysis.sh "Rachele/V2/No8-1/Epigenetic" "Design-No8-1.target" "Raw/Rachele/V2" \
                     "" "" "" \
                     "Control" "" 1.5 0.05 FALSE 1.5 0.05 FALSE v2                          
                     
bash Scripts/FullAnalysis.sh "Rachele/V2/No8-2/Epigenetic" "Design-No8-2.target" "Raw/Rachele/V2" \
                     "" "" "" \
                     "Control" "" 1.5 0.05 FALSE 1.5 0.05 FALSE v2                          
                     
bash Scripts/FullAnalysis.sh "Rachele/V1-no8/Epigenetic" "Design-no8.target" "Raw/Rachele/V1" \
                     "" "" "" \
                     "Control" "" 1.5 0.05 FALSE 1.5 0.05 FALSE v1          

bash Scripts/FullAnalysis.sh "Rachele/V1-prise2/Epigenetic" "Design-prise2.target" "Raw/Rachele/V1" \
                     "" "" "" \
                     "Control" "" 1.5 0.05 FALSE 1.5 0.05 FALSE v1                        
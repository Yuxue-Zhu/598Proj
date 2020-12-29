from timeit import default_timer as timer
from random import random
import pandas as pd
import numpy as np
from math import *
from collections import defaultdict 
from timeit import default_timer as timer

genome=open("chr22.fa","r").read()
prob=open("chr22conf.txt","r").read().split()
prob=[float(i) for i in prob]




def main():
    query1='CCCACCTCTGTGTCAAACAGTGGGGTACATGCTCTTGCTTAATCCAGCTG'
    query2='CCCATTGTGTGTCATACAGTCAGGGTACATGTTGCCTCTTAATCCAGCTG'
    #test different word length with different hsp threshold
    for i in [3,4,5,6,7]:
        w=i
        w_mer=findwords(w,genome, prob,w*log(0.5))
        r=5
        l=5
        penalty=-1
        mismatch=-1
        for z in [0.6,0.7,0.8,0.9]:
            hspthre=z
            print("w: %d, r: %d, l: %d, penalty: %d, hspthre: %f"%(w,r,l,penalty,hspthre))
            start = timer()
            alignment(query2)
            end = timer()
            print("Time: %d"%(end-start))
    #different word length with different gap penalty
    for i in [3,4,5,6,7]:
        w=i
        w_mer=findwords(w,genome, prob,w*log(0.5))
        r=5
        l=5
        mismatch=-1
        hspthre=0.7
        for z in [0.4,0.6,0.7,0.8,0.9,1,1.5,2]:
            penalty=-z
            print("w: %d, r: %d, l: %d, penalty: %f, hspthre: %f"%(w,r,l,penalty,hspthre))
            start = timer()
            alignment(query2)
            end = timer()
            print("Time: %d"%(end-start))

def alignment(query):
    hits=findhits(w_mer,query,w)
    hspScore,hsp=findHSP(w,l,r,hits,query,genome, w_mer,prob, (w+l+r)*log(hspthre))
    gapExt=getSequence(hsp,query,genome,prob,penalty, mismatch)
    q,g=getAlign(gapExt,hspScore,hsp,query,genome)

    print("acc: %f" % accur(q,g))

def accur(q,g):
    same=0
    
    for i in range(len(q)):
        
        if q[i]==g[i]:
            same+=1
    
    return same/len(q)

def getAlign(gapExt,hspScore,hsp,query,genome):
    score=float("-inf")
    hit=0
    qbest=""
    gbest=""
    ll=0
    rl=0
    for i in gapExt:

        s=0
        #print(gapExt[i])
        if gapExt[i].get("l"):
            s+=gapExt[i]["l"][2]
        s+=hspScore[i]

        if gapExt[i].get("r"):
            s+=gapExt[i]["r"][2]

        if score<s:
            score=s
            hit=i

    if gapExt[hit].get("l"):
        qbest+=gapExt[hit]["l"][0]
        gbest+=gapExt[hit]["l"][1]
        ll=gapExt[hit]["l"][3]
    qp=hsp[hit][0]
    gp=hsp[hit][1]

    qbest+=query[qp[0]:qp[1]+1]
    gbest+=genome[gp[0]:gp[1]+1]
    if gapExt[hit].get("r"):
        qbest+=gapExt[hit]["r"][0]
        gbest+=gapExt[hit]["r"][1]
        rl=gapExt[hit]["r"][3]


    print("Q: %s"%qbest)
    print("G: %s"%gbest)
    print("%d to %d" %(gp[0]-ll,gp[1]+rl))
    print("Score: %f" %score)
    return [qbest,gbest]

def NM( penalty, mismatch,genome, query, prob):

    m=[[0]*(len(query)+1) for i in range(len(genome)+1)]
    for i in range(len(genome)+1):
        m[i][0]=float('-inf') 
    for i in range(len(query)+1):
        m[0][i]=i*penalty
    for i in range(1,len(genome)+1):
        for j in range(1,len(query)+1):
            m[i][j]=max(m[i][j-1]+penalty,m[i-1][j]+penalty)
            if genome[i-1]==query[j-1]:
                m[i][j]=max(m[i][j],m[i-1][j-1]+mismatch*(1-prob[i-1]))
            else:
                m[i][j]=max(m[i][j],m[i-1][j-1]+mismatch*(2+prob[i-1])/3)
    q=""
    g=""
    i=len(genome)
    j=len(query)
    while i>0 and j>0:
        if m[i][j]-penalty==m[i][j-1]:
            g='_'+g
            q=query[j-1]+q
            j-=1
            
        elif m[i][j]-penalty==m[i-1][j]:
            g=genome[i-1]+g
            q='_'+q
            i-=1
            
        else:
       
            g=genome[i-1]+g
            q=query[j-1]+q
            i-=1
            j-=1
    if j>0:
        q=query[:j]+q
        g='_'*j+g
    return [q,g,m[len(genome)][len(query)]]

def getSequence(hsp,query,genome,prob,penalty, mismatch):
    results={}
    for s in range(len(hsp)):
        results[s]={}
        pos=hsp[s]
        qUp=query[0:pos[0][0]]
        qDn=query[pos[0][1]+1:]
        if qUp!="":
            maxs=["","",float('-inf')]
            for i in range(ceil(len(qUp)/2),2*len(qUp)+1):
                if pos[1][0]-i<0:
                    break
                else:
                    
                    upAl=NM( penalty,mismatch, genome[pos[1][0]-i:pos[1][0]], qUp, prob[pos[1][0]-i:pos[1][0]])
                    if maxs[2]<upAl[2]:
                        maxs=upAl.copy()
                        maxs.append(i)
            results[s]["l"]=maxs.copy()   
                        
        if qDn!="":
            maxs=["","",float('-inf')]
            for i in range(ceil(len(qDn)/2),2*len(qDn)+1):
                if pos[1][1]+i>=len(genome):
                    break
                else:
                    dnAl=NM( penalty,mismatch, genome[pos[1][1]+1:pos[1][1]+1+i][::-1], qDn[::-1], prob[pos[1][1]+1:pos[1][1]+1+i][::-1])
                    if maxs[2]<dnAl[2]:
                        maxs=dnAl.copy()
                        maxs[0]=maxs[0][::-1]
                        maxs[1]=maxs[1][::-1]
                        maxs.append(i)
            results[s]["r"]=maxs.copy()  
    
    return results
      

def findwords(w,genome, prob,wthre):
    words={}
    for i in range(4**w):
        word=""
        
        for j in range(w): 
            
            if i%4==0:
                word="A"+word
            elif i%4==1:
                word="T"+word
            elif i%4==2:
                word="C"+word
            else:
                word="G"+word
            i=int(i/4)
        
        words[word]={}
    for i in range(len(genome)-w+1):
 
        
        for j in list(words.keys()):
            score=0
            for z in range(w):
                if j[z]==genome[i+z]:
                    score+=logCheck(prob[i+z])
                else:
                    score+=logCheck((1-prob[i+z])/3)
    
            if score>=wthre:
                words[j][i]=score
                
    return words    


def findhits(words, query,w):
    hits=[]
    for i in range(len(query)-w+1):
        if words.get(query[i:i+w])!=[]:
            hits.append(i)
            
    return hits

def findHSP(w,l,r,hits,query,genome,w_mer, prob,ugthre):
    ungapscore=[]
    positions=[]
    for h in hits:
        
        
        for g in w_mer[query[h:h+w]]:
           # print(query[h:h+w])
            ql=h
            qr=h+w-1
            gl=g
            gr=g+w-1
            #print(ql,qr,gl,gr)
            score=w_mer[query[h:h+w]][g]
            for i in range(1,l+1):
                if h-i<0 or g-i<0:
        
                    break
                else:
                    if query[h-i]==genome[g-i]:
                        score+=logCheck(prob[g-i])
                    else:
                        score+=logCheck((1-prob[g-i])/3)
                    ql-=1
                    gl-=1
            #print(ql, gl)
            for i in range(r):
                if h+w+i>=len(query) or g+w+i>=len(genome):
                    break
                else:
                    if query[h+w+i]==genome[g+w+i]:
                        score+=logCheck(prob[g+w+i])
                    else:
                        score+=logCheck((1-prob[g+w+i])/3)
                    qr+=1
                    gr+=1
           # print(qr, gr)
          #  print(score)
            if score>=ugthre:
                positions.append([(ql,qr),(gl,gr)])
                ungapscore.append(score)
        
    return ungapscore,positions

def logCheck(i):
    if i==0:
        return float('-inf')
    else:
        return log(i)

if __name__ =='__main__':
    main()
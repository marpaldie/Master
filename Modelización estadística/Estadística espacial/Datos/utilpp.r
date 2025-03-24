# Simula n puntos aleatorios en un cuadrado unidad.
Simpalea<-function(n){matrix(runif(2*n),ncol=2)}

# Simula n quadrats circulares de radio r, aleatorios y
#       no solapados en un cuadrado unidad.
Simqalea<-function(n,r){
        if(4*r*r*n>1) stop("no caben")
        aux<-matrix(runif(2000,r,1-r),ncol=2)
        centros<-matrix(0,n,2)
        centros[1,]<-aux[1,]
        i<-1
        for(j in 2:1000){
                if(i==n) break
                pru<-aux[j,]
                if(min(sapply(1:i,function(j){dist(rbind(centros[j,],pru))}))>2*r){
                        i<-i+1
                        centros[i,]<-pru
                }
                if(j==1000) stop("demasiado complicado")
        }
        cqua<-list(centros=centros,n=n,radio=r)
        return(cqua)
}

# Dibuja los quadrats circulares dados por lqua.
dibalea.qua<-function(lqua,col="blue"){
#       points(lqua$centros,pch=21,cex=100*lqua$radio,col=col)
        theta <- seq(0, 2 * pi, length = 100)
        ro <- lqua$radio
        x <- matrix(0, nrow=lqua$n, ncol = 100)
        y <- x
        for (i in 1:lqua$n) {
            x[i, ] <- ro * cos(theta)
            y[i, ] <- ro * sin(theta)
        }
        for (i in 1:lqua$n) {
            centrosx <- rep(lqua$centros[i, 1], 100)
            centrosy <- rep(lqua$centros[i, 2], 100)
            lines(centrosx + x, centrosy + y, col = col)
        }
}

# Cuenta los puntos de un patr?n puntual pat que hay en cada uno
#       de los n quadrats circulares dados por lqua.
cupenl.qua<-function(pat,lqua){
        nump<-rep(0,lqua$n)
        for(j in 1:lqua$n)
                nump[j]<-sum(sapply(1:nrow(pat),function(i,v){dist(rbind(pat[i,],v))},
                        lqua$centros[j,])<lqua$radio)
        return(nump)
}

# Calcula varios ?ndices basados en los conteos de sucesos en quadrats
indices.qua<-function(conteos){
        suma<-sum(conteos)
        mn<-mean(conteos)
        s2<-var(conteos)
        IVR<-s2/mn
        IDM<-IVR-1
        IDL<-1+IDM/mn
        IM<-sum(conteos*(conteos-1))/(suma*(suma-1))
        cat("Indices de agrupamiento:\n")
        cat("------------------------\n")
        cat("    Varianza Relativa:     ", IVR, "\n")
        cat("    Indice de David-Moore: ", IDM, "\n")
        cat("    Desigualdad de Lloyd:  ", IDL, "\n")
        cat("    Indice de Morisita:    ", IM, "\n")
}

# Cuenta los puntos que hay en cada uno de los k*j quadrats de un grid
#       y los da en una matriz.
cuenta.qua<-function(pat,k,j){
        cortesy<-seq(0,1,len=k+1)
        cortesx<-seq(0,1,len=j+1)
        nump<-table(cut(pat[,2],cortesy),cut(pat[,1],cortesx))
        return(nump)
}

# Dibuja un grid de k*j quadrats
dibu.qua<-function(k,j,col="blue",lty="dotted"){
        cortesy<-seq(0,1,len=k+1)
        cortesx<-seq(0,1,len=j+1)
        abline(h=cortesy,col=col,lty=lty)
        abline(v=cortesx,col=col,lty=lty)
}


# Funci?n de distribuci?n de T, distancia entre dos sucesos de un patr?n
#       puntual completamente aleatorio en un cuadrado unidad.
Ht<-function(t){
        if (t<0) hh<-0
        else{
                if (t<1) hh<-pi*t*t-8*t*t*t/3+t*t*t*t/2
                else{
                        if (t<sqrt(2)) hh<-1/3-2*t*t-t*t*t*t/2+
                                4*sqrt(t*t-1)*(2*t*t+1)/3+
                                2*t*t*asin(2/(t*t)-1)
                        else hh<-1
                }
        }
        return(hh)
}

# Calcula envoltura (l?mites superior e inferior) y media de
#       simulaciones de funciones H para patrones aleatorios.
#       Hay que llamar a las librerias mva y stepfun.
#library(mva)
#library(stepfun)
Henvl<-function(nsim,pat){
        npuntos<-nrow(pat)
        ddd<-matrix(0,nsim,npuntos*(npuntos-1)/2)
        for(i in 1:nsim)
                ddd[i,]<-sort(dist(Simpalea(npuntos)))
        xu<-apply(ddd,2,min)
        xl<-apply(ddd,2,max)
        xm<-apply(ddd,2,mean)
        u<-ecdf(xu)(xu)
        l<-ecdf(xl)(xl)
        m<-ecdf(xm)(xm)
        lll<-list(xu=xu,xl=xl,xm=xm,u=u,l=l,m=m)
        return(lll)
}

# Calcula las distancias al vecino m?s pr?ximo de un patr?n puntual.
calcdvmp<-function(pat){
        ddd<-as.matrix(dist(pat))
        diag(ddd)<-rep(1000,nrow(pat))
        dvmp<-apply(ddd,1,min)
        return(dvmp)
}

# Funci?n de distribuci?n de T, distancia al vecino m?s pr?ximo de un
#       patr?n puntual completamente aleatorio en un cuadrado unidad.
Gt<-function(t,npuntos){
        if (t<0) hh<-0
        else{
                if (t<sqrt(2)) hh<-1-exp(-npuntos*pi*t*t)
                else hh<-1
        }
        return(hh)
}

# Calcula envoltura (l?mites superior e inferior) y media de
#       simulaciones de funciones G para patrones aleatorios.
#       Hay que llamar a las librerias mva y stepfun.
#library(mva)
#library(stepfun)
Genvl<-function(nsim,pat){
        npuntos<-nrow(pat)
        ddd<-matrix(0,nsim,npuntos)
        for(i in 1:nsim)
                ddd[i,]<-sort(calcdvmp(Simpalea(npuntos)))
        xu<-apply(ddd,2,min)
        xl<-apply(ddd,2,max)
        xm<-apply(ddd,2,mean)
        u<-ecdf(xu)(xu)
        l<-ecdf(xl)(xl)
        m<-ecdf(xm)(xm)
        lll<-list(xu=xu,xl=xl,xm=xm,u=u,l=l,m=m)
        return(lll)
}

# Calcula las distancias de m=k*k puntos a sus vecinos m?s
#       pr?ximos de un patr?n puntual.
caldpsmp<-function(pat,k){
        reti<-as.matrix(expand.grid(seq(1/(k+1),k/(k+1),len=k),seq(1/(k+1),k/(k+1),len=k)))
        dpsmp<-rep(0,k*k)
        for(j in 1:(k*k))
                dpsmp[j]<-min(sapply(1:nrow(pat),function(i,v){dist(rbind(pat[i,],v))},reti[j,]))
        return(dpsmp)
}

# Calcula envoltura (l?mites superior e inferior) y media de
#       simulaciones de funciones F para patrones aleatorios.
#       Hay que llamar a las librerias mva y stepfun.
#library(mva)
#library(stepfun)
Fenvl<-function(nsim,pat,k){
        npuntos<-nrow(pat)
        ddd<-matrix(0,nsim,k*k)
        for(i in 1:nsim)
                ddd[i,]<-sort(caldpsmp(Simpalea(npuntos),k))
        xu<-apply(ddd,2,min)
        xl<-apply(ddd,2,max)
        xm<-apply(ddd,2,mean)
        u<-ecdf(xu)(xu)
        l<-ecdf(xl)(xl)
        m<-ecdf(xm)(xm)
        lll<-list(xu=xu,xl=xl,xm=xm,u=u,l=l,m=m)
        return(lll)
}

myf<-function (y, trt, DFerror, MSerror, alpha = 0.05, p.adj = c("none", 
                                                            "holm", "hommel", "hochberg", "bonferroni", "BH", "BY", 
                                                            "fdr"), group = TRUE, main = NULL, console = FALSE) 
  ?match.arg('123')
match.arg(c("gauss", "rect", "ep"),
          c("gaussian", "epanechnikov", "rectangular", "triangular"),
          several.ok = TRUE)
m1<-lm(data=iris,  Sepal.Length~Species)

# 輸入 model，回傳誤差的degree of freedom
df.residual(m1)
# 輸入 model，回傳誤差的sum of square 
deviance(m1)
# MSE的算法：誤差SS/誤差自由度
MSE <-function(model) return(deviance(model)/df.residual(model))
MSE(m1)
m1 %>% anova
# pmatch: 在第二個物件裡面找尋前一個物件的元素，逐一模糊比對，並只回傳第一個找到的位置
# duplicates.ok 預設為F，意思是第二個物件被比對到之後，在下次模糊比對時，將不進行比對
pmatch("Species", names(m1$model) ) # returns 2
pmatch("med", c("mean", "median", "mode")) # returns 2
                            *
pmatch('me', c('me','me','1','me'),duplicates.ok=T) # 1

pmatch(c("", "ab", "ab"), c("abc", "ab"), duplicates.ok = TRUE)
pmatch(c("", "ab", "ab"), c("abc", "ab"), duplicates.ok = FALSE)

ipch <- pmatch(trt, names(A))


# subset 參數意思： 1:dataframe, 2:篩選條件(filter), 3.選出欄位(select)
subset(airquality, Temp > 80, select = c(Ozone, Temp))

subset(mtcars , vs==1 | am ==0 , select = c(mpg,cyl,vs,am))

# 應用：快速去掉na值data
mydata <- data.frame('x' = c(NA,NA,3,4),
                     'y' = seq(1:4) )

subset(mydata , is.na(mydata$x) == FALSE )
                       #此處

{
  p.adj <- match.arg(p.adj)
  clase <- c("aov", "lm")
  name.y <- paste(deparse(substitute(y)))
  name.t <- paste(deparse(substitute(trt)))
  if (is.null(main)) 
    main <- paste(name.y, "~", name.t)
  if ("aov" %in% class(y) | "lm" %in% class(y)) {
    if (is.null(main)) 
      main <- y$call
   
    
    # 計算MSE
    DFerror <- df.residual(y)
    MSerror <- deviance(y)/DFerror
    # A 為lm的資料表
    A <- m1$model
    
    # y 為lm的資料表中y變數的值（剛好在第一欄）
    y <- A[, 1] 
    trt = 'Species'
    names(A)[-1]
    
    # ipch 是 "處理" 在模型formula中的位置
    ipch <- pmatch(trt, names(A))
    # nipch 是 "處理" 在模型formula中的長度(有時會輸入trt為vector)
    nipch <- length(ipch)
    
    for (i in 1:nipch) {
      if (is.na(ipch[i])) 
        return(if (console) cat("Name: ", trt, "\n", 
                                names(A)[-1], "\n"))
    }
    name.t <- names(A)[ipch][1]
    
    trt <- A[, ipch]
    # 處理為1組：產生 "處理" 的data.frame
    if (nipch > 1) {
      trt <- A[, ipch[1]]
      # 處理多於兩組：產生 "處理1:處理2" 的data.frame
      for (i in 2:nipch) {
        name.t <- paste(name.t, names(A)[ipch][i], sep = ":")
        trt <- paste(trt, A[, ipch[i]], sep = ":")
      }
    }
    # y觀測值的名稱
    name.y <- names(A)[1]
  }
  ##
  junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
  ##
  Mean <- mean(junto[, 1])
  CV <- sqrt(MSerror) * 100/Mean
  medians <- tapply.stat(junto[, 1], junto[, 2], stat = "median")

  
  for (i in c(1, 5, 2:4)) {
    x <- tapply.stat(junto[, 1], junto[, 2], function(x) quantile(x)[i])
    medians <- cbind(medians, x[, 2])
  }
  
  
  medians <- medians[, 3:7]
  
  names(medians) <- c("Min", "Max", "Q25", "Q50", "Q75")
  
  means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
  sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
  nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
  std.err <- sqrt(MSerror)/sqrt(nn[, 2])
  Tprob <- qt(1 - alpha/2, DFerror)
  LCL <- means[, 2] - Tprob * std.err
  UCL <- means[, 2] + Tprob * std.err
  means <- data.frame(means, std = sds[, 2], r = nn[, 2], 
                      LCL, UCL, medians)
  names(means)[1:2] <- c(name.t, name.y)
  ntr <- nrow(means)
  nk <- choose(ntr, 2)
  
  if (p.adj != "none") {
    a <- 1e-06
    b <- 1
    for (i in 1:100) {
      x <- (b + a)/2
      xr <- rep(x, nk)
      d <- p.adjust(xr, p.adj)[1] - alpha
      ar <- rep(a, nk)
      fa <- p.adjust(ar, p.adj)[1] - alpha
      if (d * fa < 0) 
        b <- x
      if (d * fa > 0) 
        a <- x
    }
    Tprob <- qt(1 - x/2, DFerror)
  }
  
  nr <- unique(nn[, 2])
  if (console) {
    cat("\nStudy:", main)
    if (console) 
      cat("\n\nLSD t Test for", name.y, "\n")
    if (p.adj != "none") 
      cat("P value adjustment method:", p.adj, "\n")
    cat("\nMean Square Error: ", MSerror, "\n\n")
    cat(paste(name.t, ",", sep = ""), " means and individual (", 
        (1 - alpha) * 100, "%) CI\n\n")
    print(data.frame(row.names = means[, 1], means[, 2:8]))
    cat("\nAlpha:", alpha, "; DF Error:", DFerror)
    cat("\nCritical Value of t:", Tprob, "\n")
  }
  statistics <- data.frame(MSerror = MSerror, Df = DFerror, 
                           Mean = Mean, CV = CV)
  if (length(nr) == 1) 
    LSD <- Tprob * sqrt(2 * MSerror/nr)
  if (group & length(nr) == 1 & console) {
    if (p.adj == "none") 
      cat("\nleast Significant Difference:", LSD, "\n")
    else cat("\nMinimum Significant Difference:", LSD, "\n")
  }
  if (group & length(nr) != 1 & console) 
    cat("\nGroups according to probability of means differences and alpha level(", 
        alpha, ")\n")
  if (length(nr) == 1 & p.adj == "none") 
    statistics <- data.frame(statistics, t.value = Tprob, 
                             LSD = LSD)
  if (length(nr) == 1 & p.adj != "none") 
    statistics <- data.frame(statistics, t.value = Tprob, 
                             MSD = LSD)
  LSD = " "
  comb <- utils::combn(ntr, 2)
  nn <- ncol(comb)
  dif <- rep(0, nn)
  pvalue <- dif
  sdtdif <- dif
  sig <- rep(" ", nn)
  for (k in 1:nn) {
    i <- comb[1, k]
    j <- comb[2, k]
    dif[k] <- means[i, 2] - means[j, 2]
    sdtdif[k] <- sqrt(MSerror * (1/means[i, 4] + 1/means[j, 
                                                         4]))
    pvalue[k] <- 2 * (1 - pt(abs(dif[k])/sdtdif[k], DFerror))
  }
  if (p.adj != "none") 
    pvalue <- p.adjust(pvalue, p.adj)
  pvalue <- round(pvalue, 4)
  for (k in 1:nn) {
    if (pvalue[k] <= 0.001) 
      sig[k] <- "***"
    else if (pvalue[k] <= 0.01) 
      sig[k] <- "**"
    else if (pvalue[k] <= 0.05) 
      sig[k] <- "*"
    else if (pvalue[k] <= 0.1) 
      sig[k] <- "."
  }
  tr.i <- means[comb[1, ], 1]
  tr.j <- means[comb[2, ], 1]
  LCL <- dif - Tprob * sdtdif
  UCL <- dif + Tprob * sdtdif
  comparison <- data.frame(difference = dif, pvalue = pvalue, 
                           signif. = sig, LCL, UCL)
  if (p.adj != "bonferroni" & p.adj != "none") {
    comparison <- comparison[, 1:3]
  }
  rownames(comparison) <- paste(tr.i, tr.j, sep = " - ")
  if (!group) {
    if (console) {
      cat("\nComparison between treatments means\n\n")
      print(comparison)
    }
    groups <- NULL
  }
  if (group) {
    comparison = NULL
    Q <- matrix(1, ncol = ntr, nrow = ntr)
    p <- pvalue
    k <- 0
    for (i in 1:(ntr - 1)) {
      for (j in (i + 1):ntr) {
        k <- k + 1
        Q[i, j] <- p[k]
        Q[j, i] <- p[k]
      }
    }
    groups <- orderPvalue(means[, 1], means[, 2], alpha, 
                          Q, console)
    names(groups)[1] <- name.y
    if (console) {
      cat("\nTreatments with the same letter are not significantly different.\n\n")
      print(groups)
    }
  }
  parameters <- data.frame(test = "Fisher-LSD", p.ajusted = p.adj, 
                           name.t = name.t, ntr = ntr, alpha = alpha)
  rownames(parameters) <- " "
  rownames(statistics) <- " "
  rownames(means) <- means[, 1]
  means <- means[, -1]
  output <- list(statistics = statistics, parameters = parameters, 
                 means = means, comparison = comparison, groups = groups)
  class(output) <- "group"
  invisible(output)
}


library(agricolae)
treatments <- c("A","B","C")
means<-c(2,5,3)
alpha <- 0.05
pvalue<-matrix(1,nrow=3,ncol=3)
pvalue[1,2]<-pvalue[2,1]<-0.03
pvalue[1,3]<-pvalue[3,1]<-0.10
pvalue[2,3]<-pvalue[3,2]<-0.06

out<-orderPvalue(treatments,means,alpha,pvalue,console=TRUE)

barplot(out[,1],names.arg = row.names(out),col=colors()[84:87])
legend("topright",as.character(out$groups),pch=15,col=colors()[84:87],box.col=0)

install.packages('multcomp')
library(multcomp)










treatment<-treatments
n <- length(means)
z <- data.frame(treatment, means)

# 比較組字母
letras <- c(letters[1:26], LETTERS[1:26], 1:9, c(".", "+", 
                                                 "-", "*", "/", "#", "$", "%", "&", "^", "[", "]", ":", 
                                                 "@", ";", "_", "?", "!", "=", "#", rep(" ", 2000)))
# z是處理-平均值表，按平均值大到小排列
w <- z[order(z[, 2], decreasing = TRUE), ]

# M: 產生長度為處理的長度，值為''的向量
M <- rep("", n)

k <- 1
k1 <- 0
j <- 1
i <- 1

# 處理長度
# cambio:變化
cambio <- n
cambio1 <- 0
# chequeo:檢查
chequeo = 0

# 第一組為
M[1] <- letras[k]

q <- as.numeric(rownames(w))

# 分組比較的if迴圈
z

# 先判斷現在是第幾組比對，如果比對到最後一組，就停止迴圈
# 如果可以執行比對，就加1跳到下一組

while ( )

while (j < n) {
  chequeo <- chequeo + 1
  
  if (chequeo > n) 
    break
  
  for (i in j:n) {
    s <- pvalue[q[i], q[j]] > alpha
    
    if (s) {
      if (lastC(M[i]) != letras[k]) 
        M[i] <- paste(M[i], letras[k], sep = "")
    }
    
    else {
      k <- k + 1
      cambio <- i
      cambio1 <- 0
      ja <- j
      
      for (jj in cambio:n) M[jj] <- paste(M[jj], "", 
                                          sep = "")
      M[cambio] <- paste(M[cambio], letras[k], sep = "")
      
      for (v in ja:cambio) {
        if (pvalue[q[v], q[cambio]] <= alpha) {
          j <- j + 1
          cambio1 <- 1
        }
        else break
      }
      break
    }
  }
  if (cambio1 == 0) 
    j <- j + 1
}

w <- data.frame(w, stat = M)
trt <- as.character(w$treatment)
means <- as.numeric(w$means)
output <- data.frame(means, groups = M)
rownames(output) <- trt
if (k > 81) 
  cat("\n", k, "groups are estimated.The number of groups exceeded the maximum of 81 labels. change to group=FALSE.\n")
invisible(output)
}





yml() %>% yml_toc( toc_title = 'title')


# Target data is:

# Subject:
# subject | contrast | mean | ci_l | ci_h | p

# Group:
# contrast | mean | ci_l | ci_h | t | df | p

library(lmerTest)

contrast_translate <- list(
    'tot' = 'Time-on-task',
    'hr' = 'Heart rate',
    'hrv' = 'Heart rate variability',
    'power' = 'Alpha power (V^2 / Hz)',
    'fut' = 'Future thinking',
    'delib' = 'Deliberate thinking',
    'self' = 'Self-related thinking',
    'eng' = 'Disengage from thinking'
)

format_group <- function(ds, posteriors) {
    
    out <- data.frame()
    
    for (contrast in names(posteriors)) {
        print(contrast)
        post <- posteriors[[contrast]][['posterior']]$b_arou_n
        m = mean(post)
        cis <- quantile(post, probs=c(.025, .975))
        bf01 <- 1 / posteriors[[contrast]][['bf10']]
        
        # Get freq stats
        f <- as.formula(glue('{contrast} ~ arou_n + (1 + arou_n | subject)'))
        colnames(ds)[colnames(ds) == 'trial'] <- 'tot'
        if (contrast == 'power') { 
            ds <- ds[!is.na(ds$power),]
            ds$power <- log(ds$power)
        }
        mf <- lmer(f, data=ds)
        coefs <- summary(mf)$coefficients[2,]
        t <- coefs[4]
        df <- coefs[3]
        p <- coefs[5]
        sig <- ifelse(p < .05, '*', '')
        
        row <- data.frame(contrast=contrast_translate[[contrast]], 
                          mean=m, ci_l=cis[1], ci_h=cis[2], bf01=bf01,
                          t=t, df=df, p=p, sig=sig)
        out <- rbind(out, row)
    }
    
    return(out)
}

format_subjects <- function(ds, posteriors) {
    
    out <- data.frame()
    colnames(ds)[colnames(ds) == 'trial'] <- 'tot'
    
    for (contrast in names(posteriors)) {
        post <- posteriors[[contrast]][['posterior']]
        post <- post[, colnames(post) != 'b_arou_n']
        
        means <- apply(post, MARGIN=2, FUN=function (x) mean(x))
        ci_l <- apply(post, MARGIN=2, FUN=function(x) quantile(x, probs=.025))
        ci_h <- apply(post, MARGIN=2, FUN=function(x) quantile(x, probs=.975))
        subjects <- as.integer(str_extract(names(means), 'r_subject\\[(\\d+).*', group=1))
        row <- data.frame(subject=subjects, contrast=contrast_translate[[contrast]], 
                          mean=means, ci_l=ci_l, ci_h=ci_h)
        rownames(row) <- 1:nrow(row)
        out <- rbind(out, row)
    }
    
    return(out)
}

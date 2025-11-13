# Run hierarchical Bayesian models for each panel
# Output a list with names as metrics and values as sampled objects
# Load model from file if it exists

library(brms)
library(fs)
library(glue)
options(mc.cores = parallel::detectCores() - 1)
set.seed(1234)

get_samples <- function(d, contrast) {
    f_alt <- as.formula(glue('{contrast} ~ arou_n + (1 + arou_n | subject)'))
    f_null <- as.formula(glue('{contrast} ~ 1 + (1 + arou_n | subject)'))
    cores <- parallel::detectCores() - 1
    
    alt <- brm(f_alt, data=d, sample_prior='yes', 
               save_pars = save_pars(all=TRUE), iter = 5000, 
               chains=4, cores=cores)
    alt_post <- select(as_draws_df(alt), b_arou_n, matches('r_subject\\[\\d+,arou_n\\]'))
    null <- brm(f_null, data=d, sample_prior='yes', 
                save_pars = save_pars(all=TRUE), iter=5000,
                chains=4, cores=cores)
   
    rhats <- rhat(alt)
    rhat_group <- rhats['b_arou_n']
    rhat_subjects <- rhats[grepl('r_subject\\[\\d+,arou_n\\]', names(rhats))]
    
    bf10 <- bayes_factor(alt, null)
    out <- list(posterior=alt_post, bf10=bf10$bf, 
                rhat_group=rhat_group, rhat_subject=rhat_subjects)
    
    return(out)
}


get_posteriors <- function(d, overwrite=FALSE) {
    # Main function called from plot script
    
    root <- path('writing/figures/revisions/main_text/figure2_revised')
    file_path <- path(root, 'posteriors.rds')
    contrasts <- c('tot', 'hr', 'hrv', 'power', 'fut', 'delib', 'self', 'eng')
    
    # Load posteriors if they exist
    if (file.exists(file_path)) {
        result <- readRDS(file_path)
    } else result = list()
    
    
    count <- 0
    for (contrast in contrasts) {
        count <- count + 1
        print(glue('Processing contrast {count} of {length(contrasts)}'))
        # Skip already processed
        if (contrast %in% names(result) & !overwrite) {
            next
        }
        
        # Clear out missing power values (bad epochs)
        ds <- d
        colnames(ds)[colnames(ds) == 'trial'] <- 'tot'
        if (contrast == 'power') {
            ds <- ds[!is.na(ds$power),]
            # Log transform
            ds$power <- log(ds$power)
        }
        
        result[[contrast]] <- get_samples(ds, contrast)
        
        # Write out
        saveRDS(result, file_path)
        
    }
    
    return(result)
    
}

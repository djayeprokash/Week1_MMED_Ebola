# Install and load necessary packages
install.packages("remotes")
remotes::install_github("grantbrown/ABSEIR")
install.packages(c("deSolve", "ggplot2", "minpack.lm", "coda"))


library(ABSEIR)
library(deSolve)
library(ggplot2)
library(minpack.lm)
library(coda)
library(dplyr)

# Load the data
data(Kikwit1995)
head(Kikwit1995)

dim(Kikwit1995)

# Preprocess the data - Convert Date column to Date type
Kikwit1995$Date <- as.Date(Kikwit1995$Date)
# Extract the count of cases
cases <- Kikwit1995$Count
dates <- Kikwit1995$Date
times <- 1:length(cases) # Time in days

# Plot the real data
barplot(Kikwit1995$Count, xlab = "Time of Onset(Days)", ylab = "Number of New Cases",
        main = "Congo 1995")
axis(side = 1, at = seq(0, nrow(Kikwit1995), 50))

# Gaussian fitting to the real data using nlsLM
fit_gaussian <- function(data) {
  nlsLM(Count ~ a * exp(-((time - b)^2) / (2 * c^2)),
        data = data,
        start = list(a = max(data$Count), b = which.max(data$Count), c = sd(data$time)),
        control = nls.lm.control(maxiter = 100))
}

# Fit the Gaussian model
gaussian_fit <- fit_gaussian(data.frame(time = times, Count = cases))
fitted_values <- predict(gaussian_fit)
Kikwit1995$fitted <- fitted_values

# Plot the real data and Gaussian fit
ggplot() +
  geom_bar(data = Kikwit1995, aes(x = times, y = Count), stat = "identity", fill = "blue", alpha = 0.5) +
  geom_line(data = Kikwit1995, aes(x = times, y = fitted), color = "red") +
  labs(title = "Gaussian Fit to Real Data",
       x = "Time of Onset (days)", y = "Number of New Cases") +
  theme_minimal()

# Convert 'Date' to the proper format and ensure it's sorted
Kikwit1995$Date <- as.Date(Kikwit1995$Date, format = "%Y-%m-%d")
Kikwit1995 <- arrange(Kikwit1995, Date)

# Calculate the cumulative cases
Kikwit1995 <- mutate(Kikwit1995, CumulativeCases = cumsum(Count))

# Calculate days since the start of the dataset
Kikwit1995 <- mutate(Kikwit1995, DaysSinceStart = as.numeric(Date - min(Date)))

# Find the values for annotations (assuming 'Date' column exists)
mar_13_date <- as.numeric(as.Date("1995-03-13") - min(Kikwit1995$Date))
mar_13_cases <- Kikwit1995$CumulativeCases[which(Kikwit1995$DaysSinceStart == mar_13_date)]
jul_12_date <- as.numeric(as.Date("1995-07-12") - min(Kikwit1995$Date))
jul_12_cases <- Kikwit1995$CumulativeCases[which(Kikwit1995$DaysSinceStart == jul_12_date)]

ggplot(Kikwit1995, aes(x = DaysSinceStart, y = CumulativeCases)) +
   geom_line() +  # Add the line for cumulative cases
  geom_point(color="blue",shape = 1) +  # Add points for each day
  scale_x_continuous(breaks = c(0, 50, 100, 150), labels = c("0", "50", "100", "150")) +
  annotate("text", x = mar_13_date, y = mar_13_cases, label = "Mar 13", color = "white", hjust = 1.5, vjust = -0.5, fontface = "bold", size = 5) +
  annotate("text", x = jul_12_date, y = jul_12_cases, label = "Jul 12", color = "white", hjust = -0.5, vjust = -0.5, fontface = "bold", size = 5) +
  labs(title = "Epidemic Progression in Congo 1995", subtitle = "Cumulative cases plotted over days since outbreak start", x = "Days since start", y = "Cumulative number of cases") +
  theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "white"))  # Set background to white



# Define the function to calculate beta(t)
beta_t <- function(t, tau, beta_0, beta_1, q) {
  if (t < tau) {
    return(beta_0)
  } else {
    return(beta_1 + (beta_0 - beta_1) * exp(-q * (t - tau)))
  }
}


# SEIR model differential equations
seir <- function(t, y, params) {
  S <- y[1]
  E <- y[2]
  I <- y[3]
  R <- y[4]
  C <- y[5] # Cumulative cases

  beta <- beta_t(t, params["tau"], params["beta_0"], params["beta_1"], params["q"])
  N <- params["N"]
  gamma <- params["gamma"]
  k <- params["k"]

  dSdt <- -(beta * S * I) / N
  dEdt <- ((beta * S * I) / N) - (k * E)
  dIdt <- (k * E) - (gamma * I)
  dRdt <- gamma * I
  dCdt <- k * E

  return(list(c(dSdt, dEdt, dIdt, dRdt, dCdt)))
}

# Parameters and initial conditions
params <- c(
  beta_0 = 0.33,
  beta_1 = 0.09,
  k = 1/5.3,
  gamma = 1/5.61,
  tau = 56,
  q = 25,
  N = 100
)

initial_conditions <- c(S = 70, E = 29, I = 1, R = 0, C = 0)
times <- seq(1, 197, by = 1)


# Solve the SEIR model
seir_output <- ode(y = initial_conditions, times = times, func = seir, parms = params)


# Convert the output to a data frame
seir_df <- as.data.frame(seir_output)

head(seir_df)

dim(seir_df)

# Plot the simulation results
ggplot(seir_df, aes(x = time)) +
  geom_line(aes(y = I, color = "Infectious")) +
  geom_line(aes(y = E, color = "Exposed")) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = R, color = "Removed")) +
  labs(title = "SEIR Model Simulation of Ebola Outbreak",
       x = "Time (days)", y = "Number of Individuals") +
  scale_color_manual(values = c("Infectious" = "red", "Exposed" = "orange", "Susceptible" = "blue", "Removed" = "green")) +
  theme_minimal()

# Plot the real data and Gaussian fit
ggplot() +
  geom_bar(data = Kikwit1995, aes(x = times, y = Count), stat = "identity", fill = "blue", alpha = 0.5) +
  geom_line(data = seir_df, aes(x = times, y = I), color = "red") +
  labs(title = "Gaussian Fit to Real Data",
       x = "Time of Onset (days)", y = "Number of New Cases") +
  theme_minimal()

# Markov Chain SEIR model
simulate_markov_chain <- function(params, initial_conditions, times) {
  S <- initial_conditions["S"]
  E <- initial_conditions["E"]
  I <- initial_conditions["I"]
  R <- initial_conditions["R"]
  C <- initial_conditions["C"]
  
  N <- params["N"]
  beta_0 <- params["beta_0"]
  beta_1 <- params["beta_1"]
  k <- params["k"]
  gamma <- params["gamma"]
  tau <- params["tau"]
  q <- params["q"]

  results <- data.frame(time = times, S = numeric(length(times)), E = numeric(length(times)), I = numeric(length(times)), R = numeric(length(times)), C = numeric(length(times)))

  for (t in times) {
    beta <- beta_t(t, tau, beta_0, beta_1, q)
    
    # Calculate probabilities for the events
    p_SE <- 1 - exp(-(beta * I) / N)
    p_EI <- 1 - exp(-k)
    p_IR <- 1 - exp(-gamma)

    # Simulate the transitions
    new_exposed <- rbinom(1, S, p_SE)
    new_infectious <- rbinom(1, E, p_EI)
    new_removed <- rbinom(1, I, p_IR)

    # Update the compartments
    S <- S - new_exposed
    E <- E + new_exposed - new_infectious
    I <- I + new_infectious - new_removed
    R <- R + new_removed
    C <- C + new_exposed

    # Store the results
    results[results$time == t, ] <- c(t, S, E, I, R, C)
  }

  return(results)
}



# Simulate the Markov Chain model
markov_results <- simulate_markov_chain(params, initial_conditions, times)


head(markov_results)

# Plot the Markov Chain simulation results
ggplot(markov_results, aes(x = time)) +
  geom_line(aes(y = I, color = "Infectious")) +
  geom_line(aes(y = E, color = "Exposed")) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = R, color = "Removed")) +
  labs(title = "Markov Chain SEIR Model Simulation of Ebola Outbreak",
       x = "Time (days)", y = "Number of Individuals") +
  scale_color_manual(values = c("Infectious" = "red", "Exposed" = "orange", "Susceptible" = "blue", "Removed" = "green")) +
  theme_minimal()

# # Define the function to calculate R0 over time
# R0_over_time <- function(t, R0_max, beta, period) {
#   R0_max * (1 + beta * exp(-1 * pi * t / period) + cos(1.5 * pi * t / period))
# }

# # Define parameters
# R0_max <- 0.9
# beta <- 0.1
# period <- 365

# # Calculate R0 values
# R0_values <- R0_over_time(times, R0_max, beta, period)

# # Plot R0 over time
# ggplot(data.frame(time = times, R0 = R0_values), aes(x = time, y = R0)) +
#   geom_line(color = "blue") +
#   labs(title = "Reproductive Number (R0) Over Time",
#        x = "Time (days)", y = "Reproductive Number (R0)") +
#   theme_minimal()

# Calculate beta(t) dynamically based on time
beta_tt <- function(t, params) {
  if (t < params["tau"]) {
    return(params["beta_0"])
  } else {
    return(params["beta_1"] + (params["beta_0"] - params["beta_1"]) * exp(-params["q"] * (t - params["tau"])))
  }
}


# Stochastic SEIR model using Gillespie Algorithm
seir_gillespie <- function(params, initial_conditions, times) {
  t <- 0
  results <- data.frame(time = t, S = initial_conditions["S"], E = initial_conditions["E"], I = initial_conditions["I"], R = initial_conditions["R"], C = initial_conditions["C"])
  
  while(t < max(times)) {
    S <- results$S[nrow(results)]
    E <- results$E[nrow(results)]
    I <- results$I[nrow(results)]
    R <- results$R[nrow(results)]
    C <- results$C[nrow(results)]
    
    N <- params["N"]
    gamma <- params["gamma"]
    k <- params["k"]
    beta <- beta_tt(t, params)
    
    # Rates of transitions
    rate_infection <- beta * S * I / N
    rate_exposed_to_infectious <- k * E
    rate_infectious_to_removed <- gamma * I
    
    rates <- c(rate_infection, rate_exposed_to_infectious, rate_infectious_to_removed)
    total_rate <- sum(rates)
    
    if(total_rate == 0) break
    
    # Time to next event
    dt <- rexp(1, total_rate)
    t <- t + dt
    
    # Determine which event occurs
    event <- sample(1:3, 1, prob = rates)
    
    if(event == 1 && S > 0) {
      S <- S - 1
      E <- E + 1
    } else if(event == 2 && E > 0) {
      E <- E - 1
      I <- I + 1
      C <- C + 1
    } else if(event == 3 && I > 0) {
      I <- I - 1
      R <- R + 1
    }
    
    results <- rbind(results, data.frame(time = t, S = S, E = E, I = I, R = R, C = C))
  }
  
  return(results)
}


# Number of simulations
num_simulations <- 100
all_simulations <- list()

for(i in 1:num_simulations) {
  sim_result <- seir_gillespie(params, initial_conditions, times)
  all_simulations[[i]] <- sim_result
}


# Function to extract the final epidemic size from each simulation
final_epidemic_size <- function(simulation) {
  return(max(simulation$C))
}

# Calculate final epidemic sizes for all simulations
final_sizes <- sapply(all_simulations, final_epidemic_size)

# Plot the distribution of final epidemic sizes
hist(final_sizes, breaks = 30, main = "Distribution of Final Epidemic Sizes", xlab = "Final Epidemic Size", ylab = "Frequency")

# Calculate summary statistics
summary(final_sizes)


# Plot a few sample simulations
# sample_indices <- sample(1:num_simulations, 5)
sample_indices <- sample(1)


plot(NULL, xlim = c(0, max(times)), ylim = c(0, max(sapply(all_simulations, function(sim) max(sim$I)))), 
     xlab = "Time (days)", ylab = "Number of Infectious Individuals", main = "Sample Epidemic Curves")

for(i in sample_indices) {
  lines(all_simulations[[i]]$time, all_simulations[[i]]$I, col = i, type = "l")
}

legend("topright", legend = paste("Simulation", sample_indices), col = sample_indices, lty = 1)


# Parameters and initial conditions
params <- c(
  beta_0 = 0.23,
  beta_1 = 0.09,
  k = 1/5.3,
  gamma = 1/5.61,
  tau = 56,
  q = 25,
  N = 300
)

initial_conditions <- c(S = 250, E = 40, I = 10, R = 0, C = 0)
times <- seq(1, 197, by = 1)


# Stochastic SEIR model using Gillespie Algorithm
seir_gillespie <- function(params, initial_conditions, times) {
  t <- 0
  results <- data.frame(time = t, S = initial_conditions["S"], E = initial_conditions["E"], I = initial_conditions["I"], R = initial_conditions["R"], C = initial_conditions["C"])
  
  while(t < max(times)) {
    S <- results$S[nrow(results)]
    E <- results$E[nrow(results)]
    I <- results$I[nrow(results)]
    R <- results$R[nrow(results)]
    C <- results$C[nrow(results)]
    
    N <- params["N"]
    gamma <- params["gamma"]
    k <- params["k"]
    beta <- beta_tt(t, params)
    
    # Rates of transitions
    rate_infection <- beta * S * I / N
    rate_exposed_to_infectious <- k * E
    rate_infectious_to_removed <- gamma * I
    
    rates <- c(rate_infection, rate_exposed_to_infectious, rate_infectious_to_removed)
    total_rate <- sum(rates)
    
    if(total_rate == 0) break
    
    # Time to next event
    dt <- rexp(1, total_rate)
    t <- t + dt
    
    # Determine which event occurs
    event <- sample(1:3, 1, prob = rates)
    
    if(event == 1 && S > 0) {
      S <- S - 1
      E <- E + 1
    } else if(event == 2 && E > 0) {
      E <- E - 1
      I <- I + 1
      C <- C + 1
    } else if(event == 3 && I > 0) {
      I <- I - 1
      R <- R + 1
    }
    
    results <- rbind(results, data.frame(time = t, S = S, E = E, I = I, R = R, C = C))
  }
  
  return(results)
}



# Number of simulations for each tau
num_simulations <- 250

# Define range of intervention times (tau)
tau_values <- seq(0, 197, by = 1)

# Collect results
final_sizes <- data.frame(tau = tau_values, mean_final_size = NA, sd_final_size = NA)

for (tau in tau_values) {
  params["tau"] <- tau
  epidemic_sizes <- numeric(num_simulations)
  
  for (i in 1:num_simulations) {
    sim_result <- seir_gillespie(params, initial_conditions, times)
    epidemic_sizes[i] <- max(sim_result$C)
  }
  
  final_sizes[final_sizes$tau == tau, "mean_final_size"] <- mean(epidemic_sizes)
  final_sizes[final_sizes$tau == tau, "sd_final_size"] <- sd(epidemic_sizes)
}


# Load necessary library
library(ggplot2)

# Plot final epidemic size vs. time to intervention
ggplot(final_sizes, aes(x = tau, y = mean_final_size)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = mean_final_size - sd_final_size, ymax = mean_final_size + sd_final_size), alpha = 0.2) +
  labs(title = "Final Epidemic Size vs. Time to Intervention",
       x = "Time to Intervention (days)",
       y = "Final Epidemic Size") +
  theme_minimal()




# A Project for the Course STAT598W "Design & Analysis of Financial Algorithms"


# 1. Contents of the Project
In this project, I built a small Option Pricing Library with C++, including several parts as follows:
(1) Pricing of European, American, Asian, Binary, Lookback, Barrier Options.
(2) Option payoff types include: Call, Put, General payoff.
(3) Pricing methods include: BS formula, Binomial Trees, Monte Carlo Simulation, numerical
PDF.
(4) Portfolio management of stocks and options.
It includes about 1200 lines of codes in total.
And each classes and head files are explained as follows:

# 2. “stock.h”
This header file includes a stock class. It uses the Brownian Motion to update the price of
stocks. By this class, we can store basic parameters of a stock through a stock class variable, and
extract these parameters by functions inside this stock class.

# 3. “optionBS.h”
This header file includes an option class. It use BS-formula to price and updata the call and put
European options. Inside this class, it stores the underlying stock and basics parameter of this
option. And we can also extract these parameters by functions inside this option class.

# 4. “portfolio.h”
This header file includes a portfolio class. We can insert stocks or options one by one inside this
portfolio, and we can also delete stocks or options from this option. In this way, we can
management our portfolio of several stocks and option. Besides, we can get the value of this
portfolio through function inside this class.

# 5. “BinomialTree.h”
This header file includes a class of pricing option with a single underlying asset by Binomial tree
methods. It includes several functions that can pricing following options:
(1) European options with call, put and general payoff;
(2) American options with call, put and general payoff;
(3) Barrier options of call or put, up or down, in or out, (eight types)

# 6. “MonteCarlo.h”
This header file includes a class of pricing option by Monte Carlo methods. It includes several
functions that can pricing following options:
(1) European options with call and put;
(2) Binary/Digit options with call and put;
(3) Lookback options with call and put;
(4) Asian options with call and put, arithmetic mean and geometric mean;

# 7. “numericalPDE.h”
This header file includes a class of pricing option by numerical PDE methods. It uses the finite
difference method of BS PDE.
It includes a class to use explicit method of original BS PDE to price European call option.
It also includes a class to solve the PDE of Heat Diffusion Equation by explicit method.



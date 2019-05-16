# volatility-surface-constructing

This is a project about constructing a volatility surface of Bitcoin options.
There are several steps:
1. Load input files, convert the option prices into implied volatilities  
2. Build a smile for each expiry  
3. Construct a volatility surface out of the smiles
4. Create a function that asks input from user the option strike and expiry date, generate the mid, bid, ask quotes in BTC, option's Delta in BTC, and option's Vega and Theta in USD. 
Note that the mid price comes from the surface interpolation. Bid and ask prices can be BS price with a volatility +- a volatility spread say 0.5%

# Suppose we are on March 7

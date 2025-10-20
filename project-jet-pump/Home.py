import streamlit as st

st.markdown("""
<h1>🛠️ Jet Pump Performance Calculator</h1>

<p>Welcome to the <strong>Jet Pump Performance Calculator</strong> — a specialized tool designed to evaluate artificial lift efficiency in vertical oil wells using jet pump systems. This calculator integrates two industry-standard models to simulate well behavior and optimize production:</p>

<ul>
  <li><strong>Hagedorn and Brown Correlation</strong> for Vertical Lift Performance (VLP):<br>
  This empirical model estimates pressure drop in vertical multiphase flow, accounting for fluid properties, tubing geometry, and flow regime transitions. It provides a reliable bottomhole-to-wellhead pressure profile for oil injection scenarios.</li>

  <li><strong>Vogel’s Inflow Performance Relationship (IPR)</strong>:<br>
  Ideal for solution-gas drive reservoirs, Vogel’s model predicts reservoir inflow as a function of bottomhole pressure. It enables accurate estimation of production potential under varying drawdown conditions.</li>
</ul>

<h2>🔍 What This Tool Does</h2>
<ul>
  <li>Calculates and plots the <strong>VLP curve</strong> using Hagedorn and Brown for oil injection</li>
  <li>Generates the <strong>IPR curve</strong> using Vogel’s model for vertical wells</li>
  <li>Identifies the <strong>operating point</strong> where inflow meets lift capacity</li>
  <li>Supports performance diagnostics and jet pump sizing decisions</li>
</ul>

<h2>⚙️ Assumptions</h2>
<ul>
  <li>Vertical well geometry</li>
  <li>Injection fluid is <strong>oil</strong></li>
  <li>Steady-state flow conditions</li>
  <li>Reservoir drive mechanism compatible with Vogel’s IPR</li>
</ul>

<p>This calculator is built for engineers, operators, and analysts seeking fast, reliable insights into jet pump behavior and well performance. Enter your well and fluid data to begin optimizing your artificial lift strategy.</p>
""", unsafe_allow_html=True)


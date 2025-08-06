import bhmie_py

res = bhmie_py.bhmie(2.0, 1.5 + 0.1j, 90)
print(res.q_sca)

mueller = bhmie_py.mueller_mie(2.0, 1.5 + 0.1j)
print(mueller.s11)

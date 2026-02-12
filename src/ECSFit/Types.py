from collections import namedtuple

'''Used to pass various values by attribute to fit_plot().'''
lines='alka sika sikb tika1 tikb mnka1 mnkb nika nikb aula aulb auma aumb'
NamedLines = namedtuple('NamedLines', lines)

''' hold ignore limits [lower, upper] for each of the given lines.'''
LineIgnores = namedtuple('LineIgnores', 'al mn si ti')

'''XSlorentz LineE'''
LineEP = namedtuple('LineEP', 'val minoff maxoff')

'''XSlorentz width'''
WidthP = namedtuple('WidthP', 'val min max')

'''XSlorentz norm'''
NormP = namedtuple('NormP', 'val')

'''XSlorentz parameters'''
XSlorentzP = namedtuple('XSlorentzP', 'LineE width norm')

'''XSlorentz parameters for various lines'''
XSlorentzPs = namedtuple('XSlorentzPs', 'aulb aula nika nikb auma aumb sika sikb alka tika1 mnka1')


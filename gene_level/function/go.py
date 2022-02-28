import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import re
import numpy as np

def plot_GO(fname, name, term):
    
    df = pd.read_csv(fname,sep='\t')
    df = df.loc[df['Adjusted P-value']<0.05]
    if term == 'MP':
        term_name = [' '.join(segment.split(' ')[:-1]) for segment in df.Term]
    else:
        term_name = [segment.split('(')[0] for segment in df.Term]
    term_name.reverse()
    count = df.Overlap.tolist()
    count.reverse()
    bonferroni = -np.log(df['Adjusted P-value'])
    bonferroni = bonferroni.tolist()
    bonferroni.reverse()
    fig = make_subplots(
        rows=1, cols=2,
        shared_xaxes=True,
        vertical_spacing=0.03,
        specs=[[{"type": "bar"},{"type": "scatter"}]]
    )

    fig.add_trace(go.Bar(
        y=term_name,
        x=count,
        name='Count',
        orientation='h'
    ))
    
    fig.add_trace(
        go.Scatter(
            x=bonferroni,
            y=term_name,
            name="-log(padj)",
            line=dict(dash='dot')
        ),
        row=1, col=2
    )


    fig.update_layout(template='simple_white',title_text=name,xaxis = dict(showgrid=False, ticks='inside',mirror=False,showline=True),
                yaxis = dict(showgrid=False, ticks='inside',mirror=False,showline=True),yaxis2=dict(side = 'right',showticklabels=False),
                font=dict(size=18),
                  width=2000,
                  height=800,
                 title_x=0.5,
                 showlegend=True)
    fig.write_image(f'figures/{name}.pdf', scale=1, engine="kaleido")

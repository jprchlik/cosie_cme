
#add tick marks to plots
def fancy_plot(ax):
    #Turn minor ticks on
    ax.minorticks_on()
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    #set the width of the ticks
    ax.tick_params(which='both',width=1)
    #set the length of the major ticks
    ax.tick_params(which='major',length=7)
    #set length of the minor ticks
    ax.tick_params(which='minor',length=3,direction='in')
    ax.tick_params(direction='in')
    return ax

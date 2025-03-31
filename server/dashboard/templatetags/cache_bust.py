

import uuid
from django import template

register = template.Library()


@register.simple_tag(name='cache_bust')
def cache_bust():
    """
    Without this tag the horizon plot will be cached
    and not reloaded automatically after a new generation.

    https://stackoverflow.com/a/45338997
    """
    return '__v__={version}'.format(version=uuid.uuid1())

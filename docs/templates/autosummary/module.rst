{{ name | escape | underline}}

Overview
--------

.. automodule:: {{ fullname }}


   {% block attributes %}
   {% if attributes %}
   {{ _('Module Attributes') }}
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   .. autosummary::
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% if functions %}
   {{ _('Functions') }}
   ^^^^^^^^^^^^^^^^^^^^

   .. autosummary::
      :nosignatures:
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   {{ _('Classes') }}
   ^^^^^^^^^^^^^^^^^^

   .. autosummary::
      :toctree:
      :nosignatures:
      :recursive:
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   {{ _('Exceptions') }}
   ^^^^^^^^^^^^^^^^^^^^^

   .. autosummary::
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block modules %}
   {% if modules %}
   {{ _('Modules') }}
   ^^^^^^^^^^^^^^^^^^

   .. autosummary::
      :toctree:
      :recursive:
   {% for item in modules %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% if functions %}
Details
-------
.. currentmodule:: {{ fullname }}

{% for item in functions %}
.. autofunction:: {{ item }}
{%- endfor %}
{% endif %}
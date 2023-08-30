{{ name | escape | underline}}

.. currentmodule:: {{ module }}

Overview
--------

.. autoclass:: {{ objname }}

   {% block methods %}

   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :nosignatures:
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block properties %}
   {% if properties %}
   .. rubric:: {{ _('Properties') }}

   {% for item in properties %}
   .. autoproperty:: {{ objname }}.{{ item }}

   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   {% for item in attributes %}
   .. autoattribute:: {{ name }}.{{ item }}

   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   .. autosummary::
      :toctree:
      :nosignatures:
      :recursive:
   {% for item in classes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% if methods %}
Details
-------
.. class:: {{ name }}
   :noindex:
   
   {% for item in methods %}
   .. automethod:: {{ name }}.{{ item }}
   {%- endfor %}
{% endif %}
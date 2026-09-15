(function () {
  "use strict";

  function isSmallScreen() {
    return window.matchMedia && window.matchMedia("(max-width: 980px)").matches;
  }

  function currentPath() {
    try {
      return new URL(window.location.href).pathname;
    } catch (e) {
      return window.location.pathname || "/";
    }
  }

  function markActivePage() {
    var path = currentPath().replace(/\/+$/, "") || "/";
    document.querySelectorAll(".topnav a[href]:not(.nav-brand)").forEach(function (link) {
      var href = "";
      try {
        href = new URL(link.getAttribute("href"), window.location.href).pathname;
      } catch (e) {
        return;
      }
      href = href.replace(/\/+$/, "") || "/";
      if (href === path) {
        link.classList.add("active");
        link.setAttribute("aria-current", "page");
      }
    });
    var activeLink = document.querySelector(".topnav a.active");
    if (!activeLink) return;
    var group = activeLink.closest(".nav-group");
    if (!group) return;
    group.classList.add("has-active");
    var trigger = group.querySelector(".nav-group-trigger");
    if (trigger) {
      trigger.classList.add("active");
    }
  }

  function setGroupOpen(group, open) {
    if (!group) return;
    group.classList.toggle("open", Boolean(open));
    var trigger = group.querySelector(".nav-group-trigger");
    if (trigger) {
      trigger.setAttribute("aria-expanded", open ? "true" : "false");
    }
  }

  function closeGroups(nav, except) {
    if (!nav) return;
    nav.querySelectorAll(".nav-group.open").forEach(function (group) {
      if (group !== except) setGroupOpen(group, false);
    });
  }

  function setMenuOpen(nav, open) {
    if (!nav) return;
    nav.classList.toggle("nav-open", Boolean(open));
    var toggle = nav.querySelector(".nav-toggle");
    if (!toggle) return;
    toggle.setAttribute("aria-expanded", open ? "true" : "false");
    toggle.setAttribute("aria-label", open ? "收起导航" : "展开导航");
    if (!open) closeGroups(nav);
  }

  function bindNavigation() {
    var nav = document.querySelector(".topnav[data-nav]");
    if (!nav) return;

    document.addEventListener("click", function (event) {
      var target = event.target;
      if (!target || typeof target.closest !== "function") return;

      var toggle = target.closest(".nav-toggle");
      if (toggle && nav.contains(toggle)) {
        setMenuOpen(nav, !nav.classList.contains("nav-open"));
        return;
      }

      var trigger = target.closest(".nav-group-trigger");
      if (trigger && nav.contains(trigger)) {
        var group = trigger.closest(".nav-group");
        var shouldOpen = !group.classList.contains("open");
        closeGroups(nav, group);
        setGroupOpen(group, shouldOpen);
        return;
      }

      if (!nav.contains(target)) {
        setMenuOpen(nav, false);
      }
    });

    document.addEventListener("keydown", function (event) {
      if (event.key !== "Escape") return;
      setMenuOpen(nav, false);
      var toggle = nav.querySelector(".nav-toggle");
      if (toggle) toggle.focus();
    });

    window.addEventListener("resize", function () {
      if (!isSmallScreen()) {
        setMenuOpen(nav, false);
      }
    }, {passive: true});
  }

  function bindPrimaryActions() {
    document.querySelectorAll("form").forEach(function (form) {
      form.addEventListener("keydown", function (event) {
        if (event.key !== "Enter" || event.isComposing) return;
        var target = event.target;
        if (!target || !/^(INPUT|SELECT)$/.test(target.tagName)) return;
        if (target.type === "checkbox" || target.type === "radio") return;
        var primary = form.querySelector(
          "button[data-primary], button[type='submit']"
        );
        if (!primary || primary.disabled) return;
        event.preventDefault();
        primary.click();
      });
    });
  }

  function addBackToTop() {
    if (document.querySelector(".back-to-top")) return;
    var button = document.createElement("button");
    button.type = "button";
    button.className = "back-to-top";
    button.setAttribute("aria-label", "返回顶部");
    button.innerHTML = "&uarr;";
    button.addEventListener("click", function () {
      window.scrollTo({top: 0, behavior: "smooth"});
    });
    document.body.appendChild(button);

    function updateVisibility() {
      var shouldShow =
        window.scrollY > 420 &&
        (document.documentElement.scrollHeight >
          window.innerHeight + 420 || isSmallScreen());
      button.classList.toggle("visible", Boolean(shouldShow));
    }

    updateVisibility();
    window.addEventListener("scroll", updateVisibility, {passive: true});
    window.addEventListener("resize", updateVisibility, {passive: true});
  }

  function updateTaskBadge() {
    var link = document.querySelector('.topnav a[href="/tasks"]');
    if (!link) return;
    var badge = link.querySelector(".nav-count");
    if (!badge) {
      badge = document.createElement("span");
      badge.className = "nav-count";
      link.appendChild(badge);
    }
    fetch("/tasks/count", {headers: {"Accept": "application/json"}})
      .then(function (response) {
        if (!response.ok) throw new Error("bad status");
        return response.json();
      })
      .then(function (data) {
        var count = Number(data && data.count);
        if (!count || count < 0) {
          badge.hidden = true;
          badge.textContent = "";
          return;
        }
        badge.textContent = String(count);
        badge.hidden = false;
      })
      .catch(function () {
        badge.hidden = true;
      });
  }

  markActivePage();
  bindNavigation();
  bindPrimaryActions();
  addBackToTop();
  updateTaskBadge();
  window.setInterval(updateTaskBadge, 4000);
})();

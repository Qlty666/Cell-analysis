(function () {
  "use strict";

  function isSmallScreen() {
    return window.matchMedia && window.matchMedia("(max-width: 760px)").matches;
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
    document.querySelectorAll(".topnav a").forEach(function (link) {
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
    if (activeLink && isSmallScreen() && activeLink.scrollIntoView) {
      activeLink.scrollIntoView({inline: "center", block: "nearest"});
    }
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
  bindPrimaryActions();
  addBackToTop();
  updateTaskBadge();
  window.setInterval(updateTaskBadge, 4000);
})();
